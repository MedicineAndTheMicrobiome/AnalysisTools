#include "creadtagger.hpp"
#include "fourplanes.hpp"
#include "dust.hpp"

#include <limits>

USING_OLIGOFAR_SCOPES;
USING_SCOPE( fourplanes );

static Uint8 __s_debug__ = getenv( "DEBUG_CReadTagger" ) ? strtol( getenv( "DEBUG_CReadTagger" ), 0, 16 ) : 0;

CReadTagger::CReadTagger() :
    m_maxAmb( 0 ),
    m_components( 1 ),
    m_qualityChannels( 0 ),
    m_clipNwindow( 0 ),
    m_clipQuality( 0 ),
    m_clipLowercase( false ),
    m_maskComplexityEarly( false ),
    m_postLowComplexity( false ),
    m_basename( "" ),
    m_minWordsCutoff( 10 ),
    m_longWordsSwitch( 200 ),
    m_shortSeqSwitch( 101 ),
    m_negligibleSeqLength( 15 ),
    m_noPostSeqLength( 25 ),
    m_chopLength( numeric_limits<int>::max() ),
    m_chopStep( 4 ),
    m_outComplFltCutoff( 2 ),
    m_shortestWindow( numeric_limits<int>::max()/2 ),
    m_matchCountLongPct( 20,60 ),
    m_matchCountShortPct( 20,80 ),
    m_matchLongestLongPct( 10,20 ),
    m_matchLongestShortPct( 10,40 ),
    m_actions( 0 )
{}

CReadTagger::~CReadTagger()
{
    Complete();
}

void CReadTagger::Complete() 
{
    ITERATE( TBitmaskAccessVector, bma, m_bmaccess )
        if( bma->second ) delete bma->first;
    m_bmaccess.clear();
    int errors = 0;
    for( unsigned i = 0; i < sizeof( m_outputFile )/sizeof( m_outputFile[0] ); ++i ) {
        if( m_outputFile[i].get() == 0 ) continue;
        m_outputFile[i]->flush();
        m_outputFile[i]->close();
        if( m_outputFile[i]->fail() ) ++errors;
        m_outputFile[i].reset(0);
    }
    if( errors )
        THROW( runtime_error, "There was a write error for " << errors << " output file" << (errors>1?"s":"") );
}

CReadTagger::CClip CReadTagger::ClipRead( const string& read, const string& qual, const CClip& clip )
{
    int b = clip.first; 
    int e = max( b, min( clip.second, (int)read.length() ) );
    if( m_clipLowercase ) {
        while( b < e && islower( read[b] ) ) ++b;
        while( b < e && islower( read[e-1] ) ) --e;
    }
    if( m_clipQuality && qual.length() == read.length() ) {
        int qb = 0;
        int qe = read.length();
        while( qb < qe && qual[qb] - 33 < m_clipQuality ) ++qb;
        while( qb < qe && qual[qe - 1] - 33 < m_clipQuality ) --qe;
        if( qb > b ) b = qb;
        if( qe < e ) e = qe;
        if( b > e ) e = b;
    }
    if( m_clipNwindow > 0 ) {
        int bb = b;
        int ee = e;
        while( bb < e && bb - b < m_clipNwindow ) { if( toupper( read[bb] ) == 'N' ) b = bb + 1; ++bb; }
        while( b < ee && e - ee < m_clipNwindow ) { if( toupper( read[ee - 1] ) == 'N' ) e = ee - 1; --ee; }
    }
    if( e < b ) e = b;
    return CClip( b, e );
}

CReadTagger::ETagValue CReadTagger::ProcessReadData( const char * id, const char * data, const CClip& clip ) 
{
    bool newRead = false;
    if( id != m_id ) {
        PurgeRead();
        m_id = id;
        m_reads.clear();
        m_quals.clear();
        m_clips.clear();
        m_postHints.clear();
        m_decision = CDecisionMaker();
        m_forcedHuman = false;
        newRead = true;
    } else if( DecideEarly() ) return GetTagValue();

    if( data == 0 ) return m_decision.GetTagValue();
    int readLen = strlen( data );
    if( readLen == 0 ) return m_decision.GetTagValue();

    // separate read and quality
    string read, qual;
    int l = strlen( data );
    if( m_qualityChannels > 0 ) {
        ASSERT( (l&1) == 0 );
        read.assign( data, l/2 );
        qual.assign( data + l/2 );
    } else read.assign( data );

    // clip read as needed
    CClip subClip = ClipRead( read, qual, clip );
    CClip bigClip = subClip;
    ASSERT( subClip.first <= subClip.second );

    if( clip.second - clip.first < l ) {
        bigClip = ClipRead( read, qual );
        ASSERT( bigClip.first <= bigClip.second );
    }

    if( __s_debug__ & 1 ) cerr << "DEBUG:0001:" << DISPLAY( id ) << DISPLAY( read.size() ) << DISPLAY( clip ) << DISPLAY( subClip ) << DISPLAY( bigClip ) << "\n";
    if( m_actions & fAction_reportALL ) {
        ofstream& o = GetOutFile( eOutput_report );
        if( newRead ) o << id;
        if( m_actions & fAction_reportInputSeq ) o << "\t" << read;
        if( m_actions & fAction_reportInputLen ) o << "\t" << read.size();
        if( m_actions & fAction_reportClippedSeq ) o << "\t" << read.substr( bigClip.GetFrom(), bigClip.GetLength() );
        if( m_actions & fAction_reportClippedLen ) o << "\t" << bigClip.GetLength();
    }

    if( m_actions & fAction_post ) {
        m_reads.push_back( read );
        m_quals.push_back( qual );
        m_clips.push_back( clip );
        m_postHints.push_back( ePost_default );
    }
    
    if( __s_debug__ & 2 ) cerr << "DEBUG:0002:" << "\torig:\t" << read << "\n";
    if( __s_debug__ & 2 ) cerr << "DEBUG:0002:" << "\tbigc:\t" << read.substr( bigClip.GetFrom(), bigClip.GetLength() ) << "\n";
    if( __s_debug__ & 2 ) cerr << "DEBUG:0002:" << "\tclip:\t" << read.substr( clip.GetFrom(), clip.GetLength() ) << "\n";
    if( __s_debug__ & 2 ) cerr << "DEBUG:0002:" << "\tsubc:\t" << read.substr( subClip.GetFrom(), subClip.GetLength() ) << "\n";
#define SLOW 1
#if SLOW
    if( __s_debug__ & 4 ) cerr << "DEBUG:0004:" << __LINE__ << DISPLAY( m_decision ) << "\n";
    // NB: if( ....DecideEarly() ) return ...; work correct only if DecideEarly is true when and only when human decision is present
    CDecisionMaker clipped;
    clipped += ProcessClippedRead( read.c_str() + subClip.GetFrom(), subClip.GetLength(), CSeqCoding::eStrand_pos );
    if( __s_debug__ & 4 ) cerr << "DEBUG:0004:" << __LINE__ << DISPLAY( clipped ) << "\n";
    if( clipped.DecideEarly() ) return m_decision += clipped; 
    clipped += ProcessClippedRead( read.c_str() + subClip.GetFrom(), subClip.GetLength(), CSeqCoding::eStrand_neg );
    if( __s_debug__ & 4 ) cerr << "DEBUG:0004:" << __LINE__ << DISPLAY( clipped ) << "\n";
    if( clipped.DecideEarly() ) return m_decision += clipped;

    if( bigClip != subClip ) {
        CDecisionMaker unclipped;
        unclipped += ProcessClippedRead( read.c_str() + bigClip.GetFrom(), bigClip.GetLength(), CSeqCoding::eStrand_pos );
        if( __s_debug__ & 4 ) cerr << "DEBUG:0004:" << __LINE__ << DISPLAY( unclipped ) << "\n";
        if( unclipped.DecideEarly() ) return m_decision += clipped + unclipped;
        unclipped += ProcessClippedRead( read.c_str() + bigClip.GetFrom(), bigClip.GetLength(), CSeqCoding::eStrand_neg );
        if( __s_debug__ & 4 ) cerr << "DEBUG:0004:" << __LINE__ << DISPLAY( unclipped ) << "\n";
        if( unclipped.DecideEarly() ) return m_decision += clipped + unclipped;

        if( clipped.GetTagValue() == eTag_uncertain ) {
            if( unclipped.GetTagValue() == eTag_uncertain ) {
                // for this component we need to post short sequences aligned to clip as well, but post unclipped
                m_postHints.back() = ePost_both;
            } else {
                // for this component we need to post clipped sequence
                m_postHints.back() = ePost_clipped;
            }
        } else {
            // for this component we need to post unclipped version
            m_postHints.back() = ePost_unclipped;
        }
        m_decision += unclipped;
        if( __s_debug__ & 4 ) cerr << "DEBUG:0004:" << __LINE__ << DISPLAY( m_decision ) << "\n";
    }
    m_decision += clipped;
    if( __s_debug__ & 0x8 ) cerr << "DEBUG:0008:" << __LINE__ << DISPLAY( m_decision ) << "\n";
#else
    CDecisionMaker clipped = ProcessClippedRead( read.c_str() + subClip.GetFrom(), subClip.GetLength(), CSeqCoding::eStrand_pos );
    m_decision = clipped;
    if( m_decision.DecideEarly() ) return m_decision.GetTagValue();
    clipped += ProcessClippedRead( read.c_str() + subClip.GetFrom(), subClip.GetLength(), CSeqCoding::eStrand_neg );
    m_decision = clipped;
    if( bigClip != subClip ) {
        if( m_decision.DecideEarly() ) return m_decision.GetTagValue();
        CDecisionMaker unclipped = ProcessClippedRead( read.c_str() + bigClip.GetFrom(), bigClip.GetLength(), CSeqCoding::eStrand_pos );
        m_decision += unclipped;
        if( m_decision.DecideEarly() ) return m_decision.GetTagValue();
        unclipped += ProcessClippedRead( read.c_str() + bigClip.GetFrom(), bigClip.GetLength(), CSeqCoding::eStrand_neg );
        m_decision = clipped + unclipped;
        if( m_decision.DecideEarly() ) return m_decision.GetTagValue();
        if( clipped.GetTagValue() == eTag_uncertain ) {
            if( unclipped.GetTagValue() == eTag_uncertain ) {
                // for this component we need to post short sequences aligned to clip as well, but post unclipped
                m_postHints.back() = ePost_both;
            } else {
                // for this component we need to post clipped sequence
                m_postHints.back() = ePost_clipped;
            }
        } else {
            // for this component we need to post unclipped version
            m_postHints.back() = ePost_unclipped;
        }
    }
#endif
    return m_decision.GetTagValue();
}

CReadTagger::CDecisionMaker CReadTagger::ProcessClippedRead( const char * iupac, int len, CSeqCoding::EStrand strand )
{
    vector<char> read( len );
    switch( strand ) {
    case CSeqCoding::eStrand_pos: 
        for( unsigned i = 0; i < read.size(); ++i ) 
            read[i] = CNcbi8naBase( CIupacnaBase( iupac[i] ) );
        break;
    case CSeqCoding::eStrand_neg: 
        for( unsigned i = 0, j = read.size()-1; i < read.size(); ++i, --j ) 
            read[i] = CNcbi8naBase( CIupacnaBase( iupac[j] ) ).Complement();
        break;
    }
    if( GetMaskComplexityEarly() ) {
        CComplexityMeasureWindow comp; 
        int lastUnmasked = 0;
        int winsz = 0;
        ITERATE( TBitmaskAccessVector, bma, m_bmaccess ) 
            winsz = max( winsz, (int)bma->first->GetWindowLength() );
        int masked = 0;
        for( int i = 0; i < (int)read.size(); ++i ) {
            comp.Add( CNcbi8naBase(read[i]) );
            while( comp.GetWindowSize() > 66 ) comp.Del();
            if( comp.GetWindowSize() == 66 && comp.Get() > m_outComplFltCutoff ) {
                ASSERT( i >= 66-1 );
                for( lastUnmasked = max( lastUnmasked, i+winsz-1-65 ); lastUnmasked <= i-winsz+1; ++lastUnmasked ) {
                    read[lastUnmasked] = 0xf;
                    ++masked;
                }
            }
        }
        if( read.size() < 66 && comp.Get() > m_outComplFltCutoff ) {
            memset( &read[0], 0x0f, read.size() );
        }
    }
    vector<bool> bitmask;
    vector<bool> goodmask;    
    if( (int)read.size() >= m_shortestWindow ) { 
        bitmask.resize( read.size() - m_shortestWindow + 1 );
        goodmask.resize( read.size() - m_shortestWindow + 1 );
    }
    ITERATE( TBitmaskAccessVector, bma, m_bmaccess ) {
        ComputeBitmask( bitmask, goodmask, read, *bma->first, strand );
    }
    int matchLongest = 0;
    int matchCount = 0;
    int wordCount = 0;
    ComputeStatistics( bitmask, goodmask, matchLongest, matchCount, wordCount, strand );
    if( __s_debug__ & 0x0010 ) cerr << "DEBUG:0x0010" << DISPLAY( len ) << DISPLAY( strand ) << DISPLAY( read.size() ) 
        << DISPLAY( bitmask.size() ) << DISPLAY( goodmask.size() ) << DISPLAY( matchCount ) << DISPLAY( wordCount ) << DISPLAY( matchLongest ) << "\n";
    return CDecisionMaker().AddArgument( Heuristic( matchLongest, matchCount, wordCount, read.size() ) );
}

void CReadTagger::ComputeBitmask( vector<bool>& bitmask, vector<bool>& goodmask, const vector<char>& read, const CBitmaskAccess& bmaccess, CSeqCoding::EStrand )
{
    int winsz = bmaccess.GetWindowLength();
    int i0 = winsz - 1;
    if( int(read.size()) <= i0 ) { return; }
    CHashGenerator hgen( winsz );
    for( unsigned i = 0, j = -i0; i < read.size(); ++i, ++j ) {
        hgen.AddBaseMask( read[i] );
        bool found = false;
        if( hgen.GetAmbiguityCount() <= GetMaxAmbiguities() ) {
            for( CHashIterator h( hgen ); h && !found; ++h ) {
                if( bmaccess.HasWindow( *h ) ) found = true;
            }
            bitmask[j] = found;
            goodmask[j] = true;
        }
    }

}

void CReadTagger::ComputeStatistics( const vector<bool>& bitmask, const vector<bool>& goodmask, int& matchLongest, int& matchCount, int& wordCount, CSeqCoding::EStrand strand )
{
    int matchLength = matchCount = matchLongest = wordCount = 0;
    ITERATE( vector<bool>, b, goodmask ) if( *b ) ++wordCount;
    ITERATE( vector<bool>, b, bitmask ) {
        if( *b ) {
            ++matchCount;
            ++matchLength;
        } else {
            if( matchLength > matchLongest ) matchLongest = matchLength;
            matchLength = 0;
        }
    }
    if( matchLength > matchLongest ) matchLongest = matchLength;
    if( m_actions & fAction_reportALL ) {
        // here report part wil be generated
        ofstream& o = GetOutFile( eOutput_report );
        if( m_actions & fAction_reportWordGraph ) {
            o << "\t";
            if( bitmask.size() )
                ITERATE( vector<bool>, g, bitmask ) o.put( (strand == CSeqCoding::eStrand_pos ? ".>" : ".<" )[*g] );
            else 
                o << "*";
        }
        if( m_actions & fAction_reportWordTotal ) o << "\t" << wordCount;
        if( m_actions & fAction_reportWordMatchCount ) o << "\t" << matchCount;
        if( m_actions & fAction_reportWordMatchLongest ) o << "\t" << matchLongest;
    }
}

void CReadTagger::PrintHeader( EOutputFile of )
{
    switch( of ) {
        case eOutput_tag: GetOutFile( of ) << "#read-id\t#tag\n"; return;
        case eOutput_report: break;
        default: return;
    }
    ofstream& o = GetOutFile( of );
    int col = 0;
    o << "#read-id:" << ++col;
    for( int c = 1; c <= m_components; ++c ) {
        if( m_actions & fAction_reportInputSeq ) o << "\t#input-seq-" << c << ":" << ++col;
        if( m_actions & fAction_reportInputLen ) o << "\t#input-len-" << c << ":" << ++col;
        if( m_actions & fAction_reportClippedSeq ) o << "\t#clipped-seq-" << c << ":" << ++col;
        if( m_actions & fAction_reportClippedLen ) o << "\t#clipped-len-" << c << ":" << ++col;
        for( int s = 0; s <= 1; ++s ) { // by strand
            char strand = "+-"[s];
            string i = "";
            if( m_actions & fAction_reportWordGraph ) o << "\t#words-graph-" << c << strand << i << ":" << ++col;
            if( m_actions & fAction_reportWordTotal ) o << "\t#words-total-" << c << strand << i << ":" << ++col;
            if( m_actions & fAction_reportWordMatchCount ) o << "\t#words-matched-" << c << strand << i << ":" << ++col;
            if( m_actions & fAction_reportWordMatchLongest ) o << "\t#words-run-" << c << strand << i << ":" << ++col;
        }
    }
    o << "\n";
}

CReadTagger::ETagValue CReadTagger::Heuristic( int matchLongest, int matchCount, int wordCount, int readLength )
{
    if( readLength < m_negligibleSeqLength ) return eTag_foreign; // nobody cares
    if( wordCount < m_minWordsCutoff ) return eTag_uncertain;
    matchCount *= 100;
    matchLongest *= 100;
    if( wordCount >= m_longWordsSwitch ) {
        if( matchCount <  wordCount*m_matchCountLongPct.first   && matchLongest <  wordCount*m_matchLongestLongPct.first   ) return eTag_foreign;
        if( matchCount >= wordCount*m_matchCountLongPct.second  && matchLongest >= wordCount*m_matchLongestLongPct.second  ) return eTag_human;
    } else {
        if( matchCount <  wordCount*m_matchCountShortPct.first  && matchLongest <  wordCount*m_matchLongestShortPct.first  ) return eTag_foreign;
        if( matchCount >= wordCount*m_matchCountShortPct.second && matchLongest >= wordCount*m_matchLongestShortPct.second ) return eTag_human;
    }
    return eTag_uncertain;
}

CReadTagger::CDecisionMaker& CReadTagger::CDecisionMaker::operator += ( const CDecisionMaker& other ) 
{
    m_humanCnt += other.m_humanCnt;
    m_uncertainCnt += other.m_uncertainCnt;
    m_foreignCnt += other.m_foreignCnt;
    return *this;
}

CReadTagger::CDecisionMaker& CReadTagger::CDecisionMaker::AddArgument( ETagValue tvalue ) 
{ 
    switch( tvalue ) {
        case eTag_human: m_humanCnt++; break;
        case eTag_uncertain: m_uncertainCnt++; break;
        case eTag_foreign: m_foreignCnt++; break;
        default: THROW( logic_error, "Unknown tag value of " << tvalue );
    }
    return *this;
}

CReadTagger::ETagValue CReadTagger::CDecisionMaker::GetTagValue() const
{
    if( m_humanCnt ) return eTag_human;
    if( m_foreignCnt && (m_uncertainCnt == 0 ) ) return eTag_foreign;
    return eTag_uncertain;
}

bool CReadTagger::DecideEarly() const
{
    if( m_actions & fAction_reportALL ) return false;
    else if( m_forcedHuman ) return true;
    return m_decision.DecideEarly();
}

CReadTagger::ETagValue CReadTagger::GetTagValue() const 
{
    if( m_forcedHuman ) return eTag_human;
    return m_decision.GetTagValue();
}

bool CReadTagger::CDecisionMaker::DecideEarly() const
{
    if( m_humanCnt > 0 ) return true;
    else return false;
}

ofstream& CReadTagger::GetOutFile( EOutputFile of )
{ 
    static const char * suffixes[] = { ".report", ".tag", ".long.fa", ".long2.fa", ".short.fa", ".short2.fa" };
    if( unsigned(of) >= eOutput_COUNT ) THROW( logic_error, "Oops... wrong file identifier of " << of );
    if( unsigned(of) >= sizeof( m_outputFile )/sizeof( m_outputFile[0] ) ) THROW( logic_error, "Oops... wrong m_outputFile[] size of " << of );
    if( unsigned(of) >= sizeof( suffixes )/sizeof( suffixes[0] ) ) THROW( logic_error, "Oops... wrong suffixes[] size of " << of );
    if( m_outputFile[of].get() == 0 ) {
        string fname = m_basename + suffixes[of];
        m_outputFile[of].reset( new ofstream( fname.c_str() ) ); 
        if( m_outputFile[of]->fail() ) 
            THROW( runtime_error, "Failed to open file '" + fname + "': " + strerror( errno ) );
        else cerr << "* Info: Created " << fname << "\n";
        PrintHeader( of );
    }
    return *m_outputFile[of]; 
}

CReadTagger::ETagValue CReadTagger::PurgeRead()
{
    if( m_id.length() == 0 ) return eTag_uncertain;
    if( m_actions & fAction_tag ) {
        ETagValue tag = GetTagValue();
        //if( tag == eTag_uncertain && (m_actions & fAction_post) ) {}
        //else {
            ofstream& o = GetOutFile( eOutput_tag );
            o << m_id << "\t" << "FUH"[tag] << "\n"; // << "120"[tag] << "\n";
            if( o.fail() ) THROW( runtime_error, "Failed to write to tag file: " << strerror(errno) );
        //}
    }
    if( (m_actions & fAction_post) && GetTagValue() == eTag_uncertain ) {
        if( m_reads.size() >= 1 ) PrintFasta( m_reads[0], m_quals[0], 0, m_clips[0], m_postHints[0] );
        if( m_reads.size() >= 2 ) PrintFasta( m_reads[1], m_quals[0], 1, m_clips[1], m_postHints[1] );
    }
    if( m_actions & fAction_reportALL ) {
        ofstream& o = GetOutFile( eOutput_report );
        o << "\n";
        if( o.fail() ) THROW( runtime_error, "Failed to write to tag file: " << strerror(errno) );
    }
    m_reads.clear();
    m_quals.clear();
    m_clips.clear();
    m_postHints.clear();
    return GetTagValue();
}

void CReadTagger::PrintFasta( const string& read, const string& qual, int component, const CClip& clip, EPostHint hint )
{
    if( __s_debug__ & 0x20 ) cerr << "DEBUG:0020:" << DISPLAY( component ) << DISPLAY( clip ) << DISPLAY( hint ) << DISPLAY( read );
    CClip quclip( 0, read.length() );
    if( hint == ePost_clipped ) {
        CClip qcclip = ClipRead( read, qual, clip );
        PrintFasta( read.substr( qcclip.first, qcclip.GetLength() ), component, "clipped" );
    } else { //if( hint == ePost_default || hint == ePost_unclipped || hint == ePost_both ) {
        quclip = ClipRead( read, qual );
        PrintFasta( read.substr( quclip.first, quclip.GetLength() ), component, "unclipped-" + string(1,"dCUB"[hint]));
        if( hint == ePost_both ) {
            CClip qcclip = ClipRead( read, qual, clip );
            if( qcclip.second - qcclip.first > m_chopLength ) {
                if( ( qcclip.first - quclip.first ) % m_chopStep ) // then we have coord systems with different offset and should chop whole clipped read
                    PrintFasta( read.substr( qcclip.first, qcclip.GetLength() ), component, "clipped" );
                else if( (qcclip.second - quclip.first) % m_chopStep ) // we still have bad offset of last 32-mer
                    PrintFasta( read.substr( qcclip.second - m_chopLength, m_chopLength ), component, "c-tail @" + ToStr(qcclip.second - m_chopLength) );
            }
        }
    }
}

//static bool __xdebug = false;

#if 0
double CReadTagger::ComputeSimplicityIupacna( const char * iupacna ) const
{
    CComplexityMeasure dust(0);// 2 );
    Uint4 word = 0;
    int noAmb = 0;
    for( const char * c = iupacna; *c; ++c ) {
        CNcbi8naBase b8( CIupacnaBase( c ) );
        if( b8.IsAmbiguous() ) {
            dust.IncWin();
        CNcbi2naBase b2( b8.GetSmallestNcbi2na() );
        word <<= 2;
        word |= b2;
        //if( c > iupacna + 2 ) 
            dust.Add( word & 0x3f );
//        if( __xdebug ) cerr << *c << DISPLAY( c - iupacna ) << DISPLAY( dust.Get() ) << DISPLAYh( (word & 0x3f) ) << DISPLAY( dust.GetWindowSize() ) << "\n";
    }
    return dust.Get();
}
#else
double CReadTagger::ComputeSimplicityIupacna( const char * iupacna ) const
{
    CComplexityMeasure dust(0);// 2 );
    Uint4 word = 0;
    for( const char * c = iupacna, * lastAmb = c-1; *c; ++c ) {
        CNcbi4naBase b8 = CNcbi4naBase( CIupacnaBase( c ) );
        if( b8.IsAmbiguous() ) {
            dust.IncWin();
            lastAmb = c;
        } else {
            CNcbi2naBase b2( b8 ); //.GetSmallestNcbi2na() );
            word <<= 2;
            word |= b2;
            if( c - lastAmb < 3 ) dust.IncWin();
            else dust.Add( word & 0x3f );
        }
//        if( __xdebug ) cerr << *c << DISPLAY( c - iupacna ) << DISPLAY( dust.Get() ) << DISPLAYh( (word & 0x3f) ) << DISPLAY( dust.GetWindowSize() ) << "\n";
    }
    return dust.Get();
}
#endif
void CReadTagger::PrintFasta( ostream& o, const string& read, const string& comment ) 
{
    if( int(read.length()) <= m_noPostSeqLength ) return;
    double simplicity = GetPostLowComplexity() ? 0 : ComputeSimplicityIupacna( read.c_str() );
    if( simplicity <= m_outComplFltCutoff ) {
        o << ">" << m_id;
        if( comment.length() ) o << " " << comment << " ds=" << simplicity;
        o << "\n" << (read.length()?read.c_str():"n") << "\n";
//        cerr << m_id << "\t\x1b[33m" << read << "\x1b[34m\t" << simplicity << "\t\x1b[35m" << comment << "\x1b[0m\n";
        if( o.fail() ) THROW( runtime_error, "Failed to write to fasta file: " << strerror(errno) );
    }
}

void CReadTagger::PrintFasta( const string& read, int component, const string& comment )
{
    ofstream& o = GetOutFile( EOutputFile( ( int(read.length()) > m_shortSeqSwitch ? eOutput_long1 : eOutput_short1 ) + component ) );
    if( ShouldChopRead( read.length() ) ) { //<= m_shortSeqSwitch && (int)read.length() > m_chopLength ) {
        for( int i = 0; i + m_chopLength < (int)read.length(); i += m_chopStep ) {
//            __xdebug = (i == 4 || i == 65);
            PrintFasta( o, read.substr( i, m_chopLength ).c_str(), comment + " @" + ToStr( i ) );
        }
        PrintFasta( o, read.substr( read.length() - m_chopLength ).c_str(), comment + " @" + ToStr( read.length() - m_chopLength ) + "T" );
    } else {
        PrintFasta( o, read.c_str(), comment + " full" );
    }
}
