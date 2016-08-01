#ifndef cseqwcbitmask__hpp
#define cseqwcbitmask__hpp

#include "cprogressindicator.hpp"
#include "cseqvecprocessor.hpp"
#include "fourplanes.hpp"
#include "bmask-tmpl.hpp"

#include <algorithm>

BEGIN_OLIGOFAR_SCOPES

template<class BmBuilder>
class CSeqWcBitmask : public CSeqVecProcessor::ICallback
{
public:
    ~CSeqWcBitmask() {}
    CSeqWcBitmask( BmBuilder * target, Uint8 pattern = 0, int maxamb = 0, bool quiet = false ) : 
        m_target( target ), m_pattern(pattern?pattern:target->GetWordMask()), m_maxAmb(maxamb), m_quiet( quiet )
        { 
            while( m_pattern && ((m_pattern & 1) == 0) ) m_pattern >>= 1;
            m_winLength = CBitHacks::LargestBitPos( m_pattern ) + 1;
            m_pattern = CBitHacks::DuplicateBits<Uint8>( m_pattern ); 
            if( CBitHacks::BitCount8( m_pattern ) > 32 )
                THROW( logic_error, "Overrepresented word counting mode is not supported for long words " );
        }

#ifndef STANDALONE
    virtual void SequenceBegin( const TSeqIds& ids, int ) { m_seqId = ids.size() ? ids.front()->AsFastaString()  : string("unknown"); }
#else
    virtual void SequenceBegin( const TSeqIds& ids, int ) { m_seqId = ids.size() ? ids.front() : string("unknown"); }
#endif
    virtual void SequenceBuffer( CSeqBuffer * ncbi8na );
    virtual void SequenceEnd() {}

    int GetWinLength() const { return m_winLength; }

    void CommitToBitmask( int maxcount );
    typedef pair<Uint4,Uint4> TWordCount;
    typedef vector<TWordCount> TWordCounts;

protected:
    BmBuilder * m_target;
    Uint8 m_pattern;
    int   m_winLength;
    int   m_maxAmb;
    string  m_seqId;
    TWordCounts m_wordCounts;
    bool m_quiet;
};

template<class BmBuilder>
inline void CSeqWcBitmask<BmBuilder>::SequenceBuffer( CSeqBuffer * ncbi8na )
{
    CProgressIndicator progress( "Listing "+ToStr( m_target->GetWordBases() )+"-mers of length "+ToStr( GetWinLength() )+" for "+m_seqId );
    fourplanes::CHashGenerator hgen( GetWinLength() );
    vector<Uint4> words;
    words.reserve( ncbi8na->GetLength() );
    for( const char * seq = ncbi8na->GetBeginPtr(); seq < ncbi8na->GetEndPtr(); ++seq ) {
        hgen.AddBaseMask( CNcbi8naBase( seq ) );
        if( hgen.GetAmbiguityCount() <= m_maxAmb ) {
            for( fourplanes::CHashIterator h( hgen ); h; ++h ) {
                Uint8 w = CBitHacks::PackWord( h->GetLo(), m_pattern );
                ASSERT( (w&0xffffffff) == w );
                //m_target->AddExactWord( w );
                words.push_back( w );
            }
        }
        if( !m_quiet ) progress.Increment();
    }
    std::sort( words.begin(), words.end() );
    TWordCounts newcounts;
    TWordCounts::iterator x = m_wordCounts.begin(), X = m_wordCounts.end();
    vector<Uint4>::const_iterator y = words.begin(), Y = words.end();
    while( x != X && y != Y ) {
        if( x->first < *y ) ++x;
        else {
            vector<Uint4>::const_iterator n = y;
            while( n != Y && *n == *y ) ++n;
            if( x->first > *y ) 
                newcounts.push_back( make_pair( *y, n-y ) );
            else
                x->second += n-y;
            y = n;
        }
    }
    while( y != Y ) {
        vector<Uint4>::const_iterator n = y;
        while( n != Y && *n == *y ) ++n;
        newcounts.push_back( make_pair( *y, n-y ) );
        y = n;
    }
    if( newcounts.size() ) {
        m_wordCounts.reserve( m_wordCounts.size() + newcounts.size() );
        copy( newcounts.begin(), newcounts.end(), back_inserter( m_wordCounts ) );
        sort( m_wordCounts.begin(), m_wordCounts.end() );
    }
    if( !m_quiet ) progress.Summary();
}

template<class BmBuilder>
inline void CSeqWcBitmask<BmBuilder>::CommitToBitmask( int maxcount )
{
    CProgressIndicator progress( "Setting bits for words counted " + ToStr( maxcount ) + " or less" );
    for( TWordCounts::const_iterator x = m_wordCounts.begin(); x != m_wordCounts.end(); ++x ) {
        if( x->second <= maxcount && x->second > 0 ) {
            m_target->AddExactWord( x->first );
            progress.Increment();
        }
    }
    progress.Summary();
}

END_OLIGOFAR_SCOPES

#endif
