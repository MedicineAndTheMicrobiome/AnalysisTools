#ifndef OLIGOFAR_CREADTAGGER__HPP
#define OLIGOFAR_CREADTAGGER__HPP

//#include "cbitmaskaccess.hpp"
#include "bmask-tmpl.hpp"
#include "cshortreader.hpp"

BEGIN_OLIGOFAR_SCOPES

class CReadTagger
{
public:
    enum ETagValue { eTag_foreign, eTag_uncertain, eTag_human };
    enum EPostHint { ePost_default = 0, ePost_clipped = 1, ePost_unclipped = 2, ePost_both = 3 /* still don't duplicate */ };
    enum EOutputFile {
        eOutput_report,
        eOutput_tag,
        eOutput_long1,
        eOutput_long2,
        eOutput_short1,
        eOutput_short2,
        eOutput_COUNT
    };
    enum EActions {
        fAction_tag = 0x01,
        fAction_post = 0x02,
        fAction_postClipped = 0x04,
        fAction_reportInputSeq = 0x10,
        fAction_reportClippedSeq = 0x20,
        fAction_reportInputLen = 0x40,
        fAction_reportClippedLen = 0x80,
        fAction_reportWordGraph = 0x100,
        fAction_reportWordTotal = 0x200,
        fAction_reportWordMatchCount = 0x400,
        fAction_reportWordMatchLongest = 0x800,
        fAction_reportNoseq = 0xfc0,
        fAction_reportALL = 0xfff0
    };
    typedef vector<pair<CBitmaskAccess*,bool> > TBitmaskAccessVector;
    typedef Uint8 TActions;
    class CClip : public pair<int,int>
    {
    public:
        CClip( int b, int e ) : pair<int,int>( b, e ) {}
        int GetFrom() const { return first; }
        int GetToOpen() const { return second; }
        int GetTo() const { return second - 1; }
        bool Empty() const { return first >= second; }
        int GetLength() const { return second - first; }
    };
    friend ostream& operator << ( ostream& out, const CClip& clip ) { return out << "[" << clip.first << ".." << clip.second << ")"; }

    CReadTagger();
    ~CReadTagger();

    void Complete();

    int GetMinWordsCutoff() const { return m_minWordsCutoff; }
    int GetLongWordsSwitch() const { return m_longWordsSwitch; }
    int GetShortSeqSwitch() const { return m_shortSeqSwitch; }
    int GetNegligibleSeqLength() const { return m_negligibleSeqLength; }
    int GetNoPostSeqLength() const { return m_noPostSeqLength; }
    int GetChopLength() const { return m_chopLength; }
    int GetChopStep() const { return m_chopStep; }
    int GetComponents() const { return m_components; }
    int GetQualityChannels() const { return m_qualityChannels; }
    int GetMaxAmbiguities() const { return m_maxAmb; }
    int GetClipNwindow() const { return m_clipNwindow; }
    int GetClipQuality() const { return m_clipQuality; }
    bool GetClipLowercase() const { return m_clipLowercase; }
    double GetOutComplexityFilterCutoff() const { return m_outComplFltCutoff; }
    const string& GetOutputBasename() const { return m_basename; }
    const TActions& GetActions() const { return m_actions; }
    bool GetMaskComplexityEarly() const { return m_maskComplexityEarly; } 
    bool GetPostLowComplexity() const { return m_postLowComplexity; }

    int& SetMinWordsCutoff() { return m_minWordsCutoff; }
    int& SetLongWordsSwitch() { return m_longWordsSwitch; }
    int& SetShortSeqSwitch() { return m_shortSeqSwitch; }
    int& SetNegligibleSeqLength() { return m_negligibleSeqLength; }
    int& SetNoPostSeqLength() { return m_noPostSeqLength; }
    int& SetChopLength() { return m_chopLength; }
    int& SetChopStep() { return m_chopStep; }
    int& SetComponents() { return m_components; }
    int& SetQualityChannels() { return m_qualityChannels; }
    int& SetMaxAmbiguities() { return m_maxAmb; }
    int& SetClipNwindow() { return m_clipNwindow; }
    int& SetClipQuality() { return m_clipQuality; }
    bool& SetClipLowercase() { return m_clipLowercase; }
    double& SetOutComplexityFilterCutoff() { return m_outComplFltCutoff; }
    string& SetOutputBasename() { return m_basename; }
    TActions& SetActions() { return m_actions; }
    bool& SetMaskComplexityEarly() { return m_maskComplexityEarly; }
    bool& SetPostLowComplexity() { return m_postLowComplexity; }

    void SetMinWordsCutoff( int val ) { m_minWordsCutoff = val; }
    void SetLongWordsSwitch( int val ) { m_longWordsSwitch = val; }
    void SetShortSeqSwitch( int val ) { m_shortSeqSwitch = val; }
    void SetNegligibleSeqLength( int val ) { m_negligibleSeqLength = val; }
    void SetNoPostSeqLength( int val ) { m_noPostSeqLength = val; }
    void SetChopLength( int val ) { m_chopLength = val; }
    void SetChopStep( int val ) { m_chopStep = val; }
    void SetComponents( int c ) { m_components = c; }
    void SetQualityChannels( int q ) { m_qualityChannels = q; }
    void SetMaxAmbiguities( int a ) { m_maxAmb = a; }
    void SetClipNwindow( int w ) { m_clipNwindow = w; }
    void SetClipQuality( int c ) { m_clipQuality = c; }
    void SetClipLowercase( bool l ) { m_clipLowercase = l; }
    void SetOutComplexityFilterCutoff( double v ) { m_outComplFltCutoff = v; }
    void SetOutputBasename( const string& o ) { m_basename = o; }
    void SetActions( const TActions& mask, bool on ) { if( on ) m_actions |= mask; else m_actions &= ~mask; }
    void SetMaskComplexityEarly( bool on ) { m_maskComplexityEarly = on; }
    void SetPostLowComplexity( bool on ) { m_postLowComplexity = on; }

    pair<int,int> GetHeurMatchCountLongPct() const { return m_matchCountLongPct; }
    pair<int,int> GetHeurMatchCountShortPct() const { return m_matchCountShortPct; }
    pair<int,int> GetHeurMatchLongestLongPct() const { return m_matchLongestLongPct; }
    pair<int,int> GetHeurMatchLongestShortPct() const { return m_matchLongestShortPct; }

    pair<int,int>& SetHeurMatchCountLongPct() { return m_matchCountLongPct; }
    pair<int,int>& SetHeurMatchCountShortPct() { return m_matchCountShortPct; }
    pair<int,int>& SetHeurMatchLongestLongPct() { return m_matchLongestLongPct; }
    pair<int,int>& SetHeurMatchLongestShortPct() { return m_matchLongestShortPct; }

    void SetHeurMatchCountLongPct( const pair<int,int>& val ) { m_matchCountLongPct = val; }
    void SetHeurMatchCountShortPct( const pair<int,int>& val ) { m_matchCountShortPct = val; }
    void SetHeurMatchLongestLongPct( const pair<int,int>& val ) { m_matchLongestLongPct = val; }
    void SetHeurMatchLongestShortPct( const pair<int,int>& val ) { m_matchLongestShortPct = val; }

    void AddBitmaskAccess( CBitmaskAccess * bmaccess, bool takeOwnership ) {
        if( bmaccess ) { 
            m_bmaccess.push_back( make_pair( bmaccess, takeOwnership ) );
            m_shortestWindow = min( m_shortestWindow, (int)bmaccess->GetWindowLength() );
        }
    }
    
    class CDecisionMaker
    {
    public:
        CDecisionMaker() : m_humanCnt(0), m_uncertainCnt(0), m_foreignCnt(0) {}
        CDecisionMaker& operator += ( const CDecisionMaker& other );
        CDecisionMaker& AddArgument( ETagValue tvalue );
        ETagValue GetTagValue() const;
        bool DecideEarly() const;
        operator ETagValue () const { return GetTagValue(); }
        friend CDecisionMaker operator + ( const CDecisionMaker& a, const CDecisionMaker& b ) { return CDecisionMaker( a ) += b; }
        friend ostream& operator << ( ostream& o, const CDecisionMaker& m ) {
            return o << "{{h:" << m.m_humanCnt << ",u:" << m.m_uncertainCnt << ",f:" << m.m_foreignCnt << ";" << "FUH"[m.GetTagValue()] << (m.DecideEarly()?"!":"") << "}}"; 
        }
    protected:
        int m_humanCnt;
        int m_uncertainCnt;
        int m_foreignCnt;
    };

    // if qual == 1 read1, read2 are concatenamed read,qual 
    ETagValue ProcessReadData( const char * id, const char * read ) { return ProcessReadData( id, read, CClip( 0, strlen( read )/(m_qualityChannels+1) ) ); } 
    ETagValue ProcessReadData( const char * id, const char * read, int cleft, int cright ) { return ProcessReadData( id, read, CClip( cleft, cright ) ); } 
    ETagValue ProcessReadData( const char * id, const char * read, const CClip& clipOpen );
    CDecisionMaker ProcessClippedRead( const char * iupac, int len, CSeqCoding::EStrand );
    CDecisionMaker ProcessClippedRead( const vector<char>& ncbi8na, const CBitmaskAccess& bmaccess, CSeqCoding::EStrand );
    ETagValue PurgeRead();
    ETagValue GetTagValue() const;
    ETagValue Heuristic( int matchLongest, int matchCount, int wordCount, int readLength );
    void ComputeBitmask( vector<bool>& bitmask, vector<bool>& goodmask, const vector<char>& ncbi8na, const CBitmaskAccess& bmaccess, CSeqCoding::EStrand );
    void ComputeStatistics( const vector<bool>& bitmask, const vector<bool>& goodmask, int& matchLongest, int& matchCount, int& wordCount, CSeqCoding::EStrand );

    // NB: there are two types of clips: SRA clip (in clip[]) and N-,Q-,case clip. 
    // Here we use SRA clip from now on, since the other is internal and may be recomputed
    bool PostFastaClipped( EPostHint hint ) const { return m_actions & fAction_postClipped; } // || hint == ePost_clipped; }
    void ForceHuman( bool val = true ) { m_forcedHuman = val; }
    bool DecideEarly() const;

    size_t GetBitmaskCount() const { return m_bmaccess.size(); }
    const CBitmaskAccess& GetBitmask( int i ) { return *m_bmaccess[i].first; }

    CClip ClipRead( const string& read, const string& qual, const CClip& clipOpen );
    CClip ClipRead( const string& read, const string& qual ) { return ClipRead( read, qual, CClip( 0, read.length() ) ); }

    ofstream& GetOutFile( EOutputFile );
    void PrintHeader( EOutputFile );
    void PrintFasta( const string& read, const string& qual, int component, const CClip& clip, EPostHint );
    void PrintFasta( const string& read, int component, const string& comment = "" );
    void PrintFasta( ostream& o, const string& read, const string& comment = "" );
    double ComputeSimplicityIupacna( const char * iupacna ) const;
    bool ShouldChopRead( int len ) const { return len <= m_shortSeqSwitch && len > m_chopLength; }
protected:
    int m_maxAmb;
    int m_components;
    int m_qualityChannels;
    int m_clipNwindow;
    int m_clipQuality;
    bool m_clipLowercase;
    bool m_maskComplexityEarly;
    bool m_postLowComplexity;
    string m_basename;

    int m_minWordsCutoff;
    int m_longWordsSwitch;
    int m_shortSeqSwitch;
    int m_negligibleSeqLength;
    int m_noPostSeqLength;
    int m_chopLength;
    int m_chopStep;
    double m_outComplFltCutoff;

    int m_shortestWindow;
    bool m_forcedHuman;
    CDecisionMaker m_decision;
    TBitmaskAccessVector m_bmaccess;

    pair<int,int> m_matchCountLongPct;
    pair<int,int> m_matchCountShortPct;
    pair<int,int> m_matchLongestLongPct;
    pair<int,int> m_matchLongestShortPct;

    TActions m_actions;

    string m_id;
    vector<string> m_reads;
    vector<string> m_quals;
    vector<CClip > m_clips;
    vector<EPostHint> m_postHints;

    auto_ptr<ofstream> m_outputFile[eOutput_COUNT];
};

END_OLIGOFAR_SCOPES

#endif
