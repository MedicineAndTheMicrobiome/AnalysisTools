#ifndef cseqtobitmask__hpp
#define cseqtobitmask__hpp

#include "cprogressindicator.hpp"
#include "cseqvecprocessor.hpp"
#include "fourplanes.hpp"
#include "bmask-tmpl.hpp"

BEGIN_OLIGOFAR_SCOPES

template<class BmBuilder>
class CSeqToBitmask : public CSeqVecProcessor::ICallback
{
public:
    ~CSeqToBitmask() {}
    CSeqToBitmask( BmBuilder * target, Uint8 pattern = 0, int maxamb = 0, bool quiet = false ) : 
        m_target( target ), m_pattern(pattern?pattern:target->GetWordMask()), m_maxAmb(maxamb), m_quiet(quiet)
        { 
            while( m_pattern && ((m_pattern & 1) == 0) ) m_pattern >>= 1;
            m_winLength = CBitHacks::LargestBitPos( m_pattern ) + 1;
            m_pattern = CBitHacks::DuplicateBits<Uint8>( m_pattern ); 
        }
    void SetQuiet( bool to ) { m_quiet = to; }

#ifndef STANDALONE
    virtual void SequenceBegin( const TSeqIds& ids, int ) { m_seqId = ids.size() ? ids.front()->AsFastaString()  : string("unknown"); }
#else
    virtual void SequenceBegin( const TSeqIds& ids, int ) { m_seqId = ids.size() ? ids.front() : string("unknown"); }
#endif
    virtual void SequenceBuffer( CSeqBuffer * ncbi8na );
    virtual void SequenceEnd() {}

    int GetWinLength() const { return m_winLength; }

protected:
    BmBuilder * m_target;
    Uint8 m_pattern;
    int   m_winLength;
    int   m_maxAmb;
    string  m_seqId;
    bool m_quiet;
};

template<class BmBuilder>
inline void CSeqToBitmask<BmBuilder>::SequenceBuffer( CSeqBuffer * ncbi8na )
{
    CProgressIndicator progress( "Listing "+ToStr( m_target->GetWordBases() )+"-mers of length "+ToStr( GetWinLength() )+" for "+m_seqId );
    fourplanes::CHashGenerator hgen( GetWinLength() );
    for( const char * seq = ncbi8na->GetBeginPtr(); seq < ncbi8na->GetEndPtr(); ++seq ) {
        hgen.AddBaseMask( CNcbi8naBase( seq ) );
        if( hgen.GetAmbiguityCount() <= m_maxAmb ) {
            for( fourplanes::CHashIterator h( hgen ); h; ++h ) {
                Uint8 w = CBitHacks::PackWord( h->GetLo(), m_pattern );
                m_target->AddExactWord( w );
            }
        }
        if( !m_quiet ) progress.Increment();
    }
    if( !m_quiet ) progress.Summary();
}

END_OLIGOFAR_SCOPES

#endif
