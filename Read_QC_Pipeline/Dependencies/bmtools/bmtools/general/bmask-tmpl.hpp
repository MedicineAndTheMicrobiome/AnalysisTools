#ifndef OLIGOFAR_BMASK_TMPL__HPP
#define OLIGOFAR_BMASK_TMPL__HPP

#include "cbithacks.hpp"
#include "iprogress.hpp"
#include "of-debug.hpp"

#include <fstream>
#include <vector>
#include <cerrno>
#include <map>

#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

BEGIN_OLIGOFAR_SCOPES

template<class TData>
class CFlatBitmaskData
{
public:
    typedef CFlatBitmaskData<TData> self_type;
    typedef TData data_type;
    enum { kUnitBytes = sizeof(TData), kUnitBits = sizeof(TData)*8 };
    CFlatBitmaskData( void * data, int bases ) : m_data( (TData*)data ), m_bases( bases ) {}
    bool HasExactWord( Uint8 word ) const { return x_HasWord( m_data, word ); }
    void AddExactWord( Uint8 word ) { x_AddWord( m_data, word ); }
    int  GetWordBases() const { return m_bases; }
    int  GetWordBits() const { return m_bases*2; }
    Uint8 GetWordBaseMask() const { return ~((~Uint8(0)) << GetWordBases()); }
    Uint8 GetWordMask() const { return ~((~Uint8(0)) << GetWordBits()); }
    Uint8 GetDataUnits() const { return x_ComputeDataSize( m_bases ); }
    Uint8 GetDataBytes() const { return x_ComputeDataSize( m_bases )*sizeof(TData); }
    Uint8 GetDataBits() const { return x_ComputeDataSize( m_bases )*kUnitBits; }
    const TData * GetData() const { return m_data; }
    static TData * Allocate( int bases ) { return x_Allocate( bases ); }
    static Uint8 GetSize( int bases ) { return x_ComputeDataSize( bases ); }
    static Uint8 GetBits( int bases ) { return Uint8(1) << (bases*2); }
    static Uint8 GetBytes( int bases ) { return GetBits( bases )/8; }
    static Uint8 GetUnitBits() { return kUnitBits; }
    class C_Iterator
    {
    public:
        C_Iterator( const CFlatBitmaskData& data, Uint8 word = 0 ) : 
            m_data( data.GetData() ), m_openLimit( data.GetWordMask() + 1 ), 
            m_word(x_GetFirstWord(data.GetData(),word,m_openLimit)) { }
        Uint8 operator * () const { return m_word; }
        const Uint8 * operator -> () const { return &m_word; }
        C_Iterator& operator ++ () { m_word = x_GetNextWord(m_data,m_word,m_openLimit); return *this; }
        C_Iterator operator ++ ( int ) { C_Iterator x( *this ); ++*this; return x; } 
        bool operator == ( const C_Iterator& other ) const { return m_word == other.m_word; }
        bool operator != ( const C_Iterator& other ) const { return m_word != other.m_word; }
    protected:
        const TData * m_data;
        Uint8 m_openLimit;
        Uint8 m_word;
    };
    typedef C_Iterator const_iterator;
    typedef Uint8 value_type;
    const_iterator begin() const { return C_Iterator( *this, 0 ); }
    const_iterator end() const { return C_Iterator( *this, GetWordMask()+1 ); }
    template<class Callback>
    void ForEachSet( Callback& cbk ) const { 
        Uint8 W = (1ULL << (m_bases*2)); 
        //DEBUG( DISPLAY( W ) << DISPLAY( *begin() ) << DISPLAY( *end() ) );
        for( const_iterator x = begin(), X = end(); x != X; ++x ) { ASSERT( *x < W ); cbk( *x ); } 
    }
protected:
    static bool x_HasWord( const TData * data, Uint8 word ) {
        return data[word / kUnitBits] & (1UL << (word % kUnitBits));
    }
    static void x_AddWord( TData * data, Uint8 word ) {
        data[word / kUnitBits] |= (1UL << (word % kUnitBits));
    }
    static Uint8 x_ComputeDataSize( int wordBases ) {
        return (1ULL << (2*wordBases))/sizeof(TData)/8;
    }
    static Uint8 x_GetFirstWord( const TData * data, Uint8 word, Uint8 openLimit ) {
        if( word >= openLimit ) return openLimit;
        if( x_HasWord(data, word) ) return word;
        else return x_GetNextWord( data, word, openLimit );
    }
    static Uint8 x_GetNextWord( const TData * data, Uint8 word, Uint8 openLimit ) {
        if( word >= openLimit ) return openLimit;
        do ++word; while( word < openLimit && !x_HasWord( data, word ) );
        return word;
    }
    static TData * x_Allocate( int wordBases ) {
        Uint8 size = x_ComputeDataSize( wordBases );
        TData * data = new TData[ size ];
        memset( data, 0, sizeof( data[0] )*size );
        //DEBUG("Allocated " << size << " units of size " << sizeof( Uint8 ) << " at " << data);
        return data;
    }
protected:
    TData * m_data;
    int     m_bases;
};

template<Uint8 value> class CUint8const { public: typedef Uint8 value_type; static Uint8 Get() { return value; } };
template<Uint4 value> class CUint4const { public: typedef Uint4 value_type; static Uint4 Get() { return value; } };
template<int value> class CIntConst { public: typedef int value_type; static int Get() { return value; } };
template<typename TValue> class CIntVal 
{ 
public: 
    typedef TValue value_type; 
    void Set( const TValue& val ) { m_data = val; }
    TValue  Get() const { return m_data; } 
    TValue& Set() { return m_data; }
    CIntVal( TValue val ) : m_data( val ) {} 
protected: 
    TValue m_data; 
};

template<class TPrefixBits, class TSuffixBits, class TOffsetBits, class TCountBits>
class CBmPackConfig
{
public:
    CBmPackConfig( const TPrefixBits& pfxBits, const TSuffixBits& sfxBits, const TOffsetBits& offBits, const TCountBits & cntBits ) :
        m_pfxBits( pfxBits ),
        m_sfxBits( sfxBits ),
        m_offBits( offBits ),
        m_cntBits( cntBits ) { ASSERT( GetCountMask() == 0 || GetLongestList() < GetCountMask() ); }
    typename TPrefixBits::value_type  GetPrefixBits() const { return m_pfxBits.Get(); }
    typename TSuffixBits::value_type  GetSuffixBits() const { return m_sfxBits.Get(); }
    typename TOffsetBits::value_type  GetOffsetBits() const { return m_offBits.Get(); }
    typename TCountBits ::value_type  GetCountBits() const { return m_cntBits.Get(); }
    typename TPrefixBits::value_type  GetWordBits() const { return m_pfxBits.Get() + m_sfxBits.Get(); }
    typename TOffsetBits::value_type  GetRecordBits() const { return GetCountBits() + GetOffsetBits(); }
    typename TOffsetBits::value_type  GetRecordBytes() const { return (GetRecordBits()+7)/8; }

    Uint8 GetRecordCount() const { return 1UL << GetPrefixBits(); }
    Uint8 GetTableBytesPacked() const { return (GetRecordCount()*GetRecordBits()+7)/8; }
    Uint8 GetTableBytesAligned() const { return GetRecordCount()*GetRecordBytes(); }
    Uint8 GetLongestList() const { return GetSuffixBits()?(1UL << GetSuffixBits())/GetSuffixBits():0; }
    Uint8 GetBitmaskBits() const { return (1UL<<GetSuffixBits()); }
    Uint8 GetBitmaskBytes() const { return GetSuffixBits()<=3?1:(1UL<<(GetSuffixBits()-3)); }
    Uint8 GetListBits( unsigned count ) const { return GetSuffixBits()*count; }
    Uint8 GetListBytes( unsigned count ) const { return (GetListBits(count)+7)/8; }
    Uint8 GetDataEntryBits( unsigned count ) const { return count > GetLongestList() ? GetBitmaskBits() : GetListBits( count ); }
    Uint8 GetDataEntryBytes( unsigned count ) const { return count > GetLongestList() ? GetBitmaskBytes() : GetListBytes( count ); }
    Uint8 GetWordBases() const { return (GetWordBits() + 1)/2; }
    Uint8 GetWordBaseMask() const { return ~((~Uint8(0))<<GetWordBases()); }
    Uint8 GetWordMask() const { return ~((~Uint8(0))<<GetWordBits()); }
    Uint8 GetIndexMask() const { return ~((~Uint8(0))<<GetPrefixBits()); }
    Uint8 GetCountMask() const { return ~((~Uint8(0))<<GetCountBits()); }
    Uint8 GetSuffixMask() const { return ~((~Uint8(0))<<GetSuffixBits()); }
    Uint8 GetPrefixMask() const { return GetIndexMask()<<GetSuffixBits(); }
    Uint8 GetOffsetMask() const { return ~((~Uint8(0))<<GetOffsetBits()); }
    Uint8 Word2Word(Uint8 word) const { return word & GetWordMask(); }
    Uint8 Word2Prefix(Uint8 word) const { return word&GetPrefixMask(); }
    Uint8 Word2Suffix(Uint8 word) const { return word&GetSuffixMask(); }
    Uint8 Word2Index(Uint8 word) const { return (word >> GetSuffixBits()); }
protected:
    TPrefixBits m_pfxBits;
    TSuffixBits m_sfxBits;
    TOffsetBits m_offBits;
    TCountBits  m_cntBits;
};


template <class TPackConfig>
class CPackedBitmaskData
{
public:
    CPackedBitmaskData( const TPackConfig& pcfg, void * table, void * data ) :
        m_packCfg( pcfg ), m_table( (Uint1*)table ), m_data( (Uint1*)data ) {}
    const TPackConfig& GetConfig() const { return m_packCfg; }
    const void * GetData()  const { return m_data; }
    const void * GetTable() const { return m_table; }
protected:
    Uint8 GetSuffixCount( Uint8 index ) const { return CBitHacks::GetBits( m_table, index*m_packCfg.GetRecordBits(), m_packCfg.GetCountBits() ); }
    Uint8 GetDataOffset( Uint8 index ) const { return  CBitHacks::GetBits( m_table, index*m_packCfg.GetRecordBits() + m_packCfg.GetCountBits(), m_packCfg.GetOffsetBits() ); }
    Uint8 GetTableEntry( Uint8 offset, int listOrd ) const { return CBitHacks::GetBits( m_data + offset, listOrd*m_packCfg.GetSuffixBits(), m_packCfg.GetSuffixBits() ); }
    bool  GetBitmaskValue( Uint8 offset, Uint8 suffix ) const { return CBitHacks::GetBits( m_data + offset, suffix, 1 ); }
protected:
    TPackConfig m_packCfg;
    Uint1 * m_table;
    Uint1 * m_data;
};

template <class TPackConfig>
class CPackedBitmaskDataReader : public CPackedBitmaskData<TPackConfig>
{
public:
    typedef Uint8 value_type;
    typedef CPackedBitmaskData<TPackConfig> super_type;
    typedef CPackedBitmaskDataReader<TPackConfig> self_type;
    CPackedBitmaskDataReader( const TPackConfig& pcfg, const void * table, const void * data ) :
        CPackedBitmaskData<TPackConfig>( pcfg, const_cast<void*>(table), const_cast<void*>(data) ) {}
    bool HasExactWord( Uint8 word ) const;
    template<class Callback>
    void ForEachSet( Callback& cbk ) const;

    template<class Iterator>
    Iterator GetSuffixList( Uint8 word, Iterator dest ) const;
    Uint8 GetSuffixListCount( Uint8 word ) { return GetSuffixCount( word >> super_type::GetConfig().GetSuffixBits() ); }
};

template<class TPackConfig>
class CPackedBitmaskDataAppender : public CPackedBitmaskData<TPackConfig>
{
public:
    typedef Uint8 value_type;
    typedef vector<Uint4> TSuffixes;
    typedef CPackedBitmaskData<TPackConfig> super_type;
    typedef CPackedBitmaskDataAppender<TPackConfig> self_type;
    typedef vector<Uint8> TOffsets;
    typedef map<Uint4,TOffsets> TOffsetHash;

    ~CPackedBitmaskDataAppender() { Free(); }
    CPackedBitmaskDataAppender( const TPackConfig& pcfg ) :
        super_type( pcfg, new Uint1[pcfg.GetTableBytesPacked()], 0 ), m_dataSize(0), m_lastWord(0), m_lastIndex(0), m_bitCount(0), m_addCount(0), m_savedBytes(0), m_useHashOffsets( false ) {}
    void Init() { Free(); super_type::m_table = new Uint1[super_type::GetConfig().GetTableBytesPacked()]; }
    void Free() { delete[] super_type::m_table; free( super_type::m_data ); Release(); }
    void Release() { super_type::m_table = 0; super_type::m_data = 0; m_dataSize = 0; m_lastIndex = 0; m_lastWord = 0; m_bitCount = m_addCount = m_savedBytes = 0; m_suffixes.clear(); m_offsetHash.clear(); }
    bool Valid() const { return super_type::m_table != 0; }
    void AddExactWord( Uint8 word );
    void Complete() { Purge(); }
    Uint8 GetDataSize() const { return m_dataSize; }
    Uint8 GetBitCount() const { return m_bitCount; }
    Uint8 GetAddCount() const { return m_addCount; }
    Uint8 GetSavedBytes() const { return m_savedBytes; }
    void SetUseHashOffsets( bool value ) { m_useHashOffsets = value; }
protected:
    void Purge();
protected:
    Uint8 m_dataSize;
    Uint8 m_lastWord;
    Uint8 m_lastIndex;
    Uint8 m_bitCount;
    Uint8 m_addCount;
    Uint8 m_savedBytes;
    TSuffixes m_suffixes;
    TOffsetHash m_offsetHash;
    bool m_useHashOffsets;
};

typedef CBmPackConfig<CIntVal<Uint2>,CIntVal<Uint2>,CIntVal<Uint2>,CIntVal<Uint2> > CVarBmPackConfig;

class CBitmaskFileHeader : public CVarBmPackConfig
{
public:
    bool IsCompressed() const { return GetPrefixBits(); }
    void Read( ifstream& );
    void Write( ofstream& );
    void Print( ostream& );
    void Reconsile( int minor );
    void SetFileMinor( int ver );

    void SetPattern( Uint8 basePattern );
    void SetPrefixBits( int bits ); // use 0 for flat file
    void SetUnitBits( int bits ) { m_bitsPerUnit = bits; m_unitsTotal = (GetWordMask()+1)/bits; }
    void SetWordsPresent( Uint8 words ) { m_wordsPresent = words; }
    void SetDataSize( Uint8 bytes ) { m_dataSize = bytes; }

    unsigned GetSignatureLength() const { return m_signatureLength; }
    unsigned GetSignatureSize() const { return m_signature.size(); }
    const char * GetSignature() const { return &m_signature[0]; }

    int GetFileMajor()   const { return m_signature[m_signatureLength+1]-'0'; }
    int GetFileMinor()   const { return m_signature[m_signatureLength+3]-'0'; }
    int GetFileRelease() const { return m_signature[m_signatureLength+5]-'0'; }

    Uint8 GetByteOrder() const { return m_byteOrder; }
    Uint4 GetHeaderSize() const { return m_headerSize; }
    Uint4 GetBitsPerUnit() const { return m_bitsPerUnit; }
    Uint4 GetUnitsTotal() const { return m_unitsTotal; }
    Uint4 GetMaxAmb() const { return m_maxAmb; }
    Uint4 GetWinStep() const { return m_winStep; }
    Uint8 GetWordsPresent() const { return m_wordsPresent; }
    Uint4 GetBasePattern() const { return m_basePattern; }
    Uint8 GetNcbi2naPattern() const { return m_ncbi2naPattern; }
    //Uint8 GetWordMask() const { return m_wordMask; }
    Uint2 GetWindowLength() const { return m_winLength; }
    Uint8 GetFlags() const { return m_flags; }
    Uint8 GetContentOffset() const { return m_contentOffset; }
    Uint8 GetTableOffset() const { return m_contentOffset; }
    Uint8 GetDataOffset() const { return m_contentOffset + GetTableBytesPacked(); }
    Uint8 GetDataSize() const { return m_dataSize; }
    CBitmaskFileHeader( const CVarBmPackConfig& cfg = CVarBmPackConfig(0,0,0,0) );

    bool IsDiscontiguous() const { return GetWindowLength() > GetWordBases(); }
    bool IsContiguous() const { return GetWindowLength() == GetWordBases(); }
protected:
    void x_InitSignature();
    template <int szof,class T>
    static void x_Read( istream& in, T* dest, int cnt = 1 );
    template <int szof,class T>
    static void x_Write( ostream& out, const T* src, int cnt = 1 );
protected:
    Uint8 m_byteOrder;
    Uint4 m_headerSize;
    Uint4 m_bitsPerUnit;
    Uint8 m_unitsTotal;
    Uint4 m_maxAmb;
    Uint4 m_winStep;
    Uint8 m_wordsPresent;
    Uint4 m_basePattern;
    Uint8 m_ncbi2naPattern;
    //Uint8 m_wordMask;
    Uint2 m_winLength;
    Uint8 m_flags;
    Uint8 m_contentOffset;
    Uint8 m_dataSize;
    vector<char> m_signature;
    int m_signatureLength;
};

class CBitmaskAccess : public CBitmaskFileHeader
{
public:
    enum EAttachFlags { fUseMmap = 1, fMmapSequential = 2 };
    typedef CFlatBitmaskData<Uint4> TFlatBitmap;
    typedef CPackedBitmaskDataReader<CVarBmPackConfig> TPackedBitmap;
    
    CBitmaskAccess();
    ~CBitmaskAccess() { ReleaseFile(); }

    void AttachFile( const string& filename, int flags, IProgressIndicator * p = 0 );
    void ReleaseFile();

    bool IsFlat() const { return !m_packedBm.GetConfig().GetPrefixBits(); }
    bool IsPacked() const { return m_packedBm.GetConfig().GetPrefixBits(); }

    const TFlatBitmap& GetFlatBitmask() const { return m_flatBm; }
    const TPackedBitmap& GetPackedBitmask() const { return m_packedBm; }
    
    bool HasExactWord( Uint8 word ) const { if( IsPacked() ) return m_packedBm.HasExactWord( word ); else return m_flatBm.HasExactWord( word ); }
    bool HasWindow( Uint8 window ) const { Uint8 word = CBitHacks::PackWord( window, GetNcbi2naPattern() ); return HasExactWord( word ); }

    template<class Target>
    class C_Builder 
    {
    public:
        C_Builder( Target& t, IProgressIndicator * p ) : m_target( t ), m_progress( p ) {}
        void operator () ( Uint8 word ) { m_target.AddExactWord( word ); if( m_progress ) m_progress->Increment(); }
    protected:
        Target& m_target;
        IProgressIndicator * m_progress;
    };
    
    template<class Callback>
    void ForEachSet( Callback& cbk ) const {
       if( IsPacked() ) m_packedBm.ForEachSet( cbk ); else m_flatBm.ForEachSet( cbk );
    }
protected:
    TFlatBitmap m_flatBm;
    TPackedBitmap m_packedBm;
    unsigned m_xflags;
};

class CBitmaskAccessExt : public CBitmaskAccess
{
public:
    void SetRqWordSize( int sz ) { m_rqWordSize = sz; }
    int  GetRqWordSize() const { return m_rqWordSize; }
    CBitmaskAccessExt() : m_rqWordSize(0) {}
    bool HasWord( Uint8 word ) const;
protected:
    int m_rqWordSize;
};

class CBitmaskWriter : public CBitmaskFileHeader
{
public:
    CBitmaskWriter() : m_maxAmb(0), m_minor(0) {}
    class CCounter
    {
    public:
        CCounter() : m_count(0) {}
        void operator () ( Uint8 ) { ++m_count; }
        Uint8 GetCount() const { return m_count; }
    protected:
        Uint8 m_count;
    };
    void SetFileMinor( int m ) { m_minor = m; }

    void SetMaxAmb( int ma ) { m_maxAmb = ma; }
    template<class Unit>
    void Write( const string& fname, const CFlatBitmaskData<Unit>& data, IProgressIndicator * p = 0 );
    template<class Config>
    void Write( const string& fname, const CPackedBitmaskData<Config>& data, Uint8 dsize, IProgressIndicator * p = 0 );
    template<class TData>
    void InitHeader( CBitmaskFileHeader& header, const TData& data ) { InitHeader( header, data, data ); }
    template<class TData, class TCfg>
    void InitHeader( CBitmaskFileHeader& header, const TData& data, const TCfg& cfg );
protected:
    int m_maxAmb;
    int m_minor;
};

////////////////////////////////////////////////////////////////////////
// Implementations

template<class TPackConfig>
inline bool CPackedBitmaskDataReader<TPackConfig>::HasExactWord( Uint8 word ) const 
{
    word = super_type::GetConfig().Word2Word( word );
    Uint8 index = super_type::GetConfig().Word2Index( word );
    Uint8 count = super_type::GetSuffixCount( index );
    if( count == 0 ) return false;
    Uint8 table = super_type::GetDataOffset( index );
    if( count > super_type::GetConfig().GetLongestList() ) {
        return CBitHacks::GetBits( super_type::m_data + table, super_type::GetConfig().Word2Suffix( word ), 1 );
    } else if( count < (Uint8)super_type::GetConfig().GetSuffixBits()/2 ) {
        Uint8 wordsfx = super_type::GetConfig().Word2Suffix( word );
        for( Uint4 x = 0; x < count; ++x ) {
            Uint8 suffix = super_type::GetTableEntry( table, x );
            if( suffix == wordsfx ) return true;
        } 
    } else {
        Uint8 wordsfx = super_type::GetConfig().Word2Suffix( word );
        Int8 a = 0, b = count - 1, x = count/2;
        while( a <= b ) {
            Uint8 suffix = super_type::GetTableEntry( table, x );
            if( suffix > wordsfx ) b = x - 1;
            else if( suffix < wordsfx ) a = x + 1;
            else return true;
            x = (a + b)/2;
        }
    }
    return false;
}

template<class TPackConfig>
template<class Iterator>
inline Iterator CPackedBitmaskDataReader<TPackConfig>::GetSuffixList( Uint8 word, Iterator dest ) const
{
    word = super_type::GetConfig().Word2Word( word );
    Uint8 index = super_type::GetConfig().Word2Index( word );
    Uint8 count = super_type::GetSuffixCount( index );
    if( count == 0 ) return dest;
    Uint8 table = super_type::GetDataOffset( index );
    if( count > super_type::GetConfig().GetLongestList() ) {
        for( Uint8 sfx = 0; sfx <= super_type::GetConfig().GetSuffixMask(); ++sfx ) {
            if( CBitHacks::GetBits( super_type::m_data + table, sfx, 1 ) ) *dest++ = sfx;
        }
    } else {
        for( Uint4 x = 0; x < count; ++x ) {
            Uint8 suffix = super_type::GetTableEntry( table, x );
            *dest++ = suffix;
        } 
    }
    return dest;
}



template<class TPackConfig>
template<class Callback>
void CPackedBitmaskDataReader<TPackConfig>::ForEachSet( Callback& cbk ) const
{
    for( Uint8 index = 0; index <= super_type::GetConfig().GetIndexMask(); ++index ) {
        Uint8 count = super_type::GetSuffixCount( index );
        if( count == 0 ) continue;
        Uint8 offset = super_type::GetDataOffset( index );
        Uint8 prefix = index << super_type::GetConfig().GetSuffixBits();
        if( count > super_type::GetConfig().GetLongestList() ) {
            for( Uint8 suffix = 0; suffix <= super_type::GetConfig().GetSuffixMask(); ++suffix ) {
                if( super_type::GetBitmaskValue( offset, suffix ) ) cbk( prefix | suffix );
            }
        } else {
            for( Uint8 sno = 0; sno < count; ++sno ) {
                cbk( prefix | super_type::GetTableEntry( offset, sno ) );
            }
        }
    }
}

template<class TPackConfig>
inline void CPackedBitmaskDataAppender<TPackConfig>::AddExactWord( Uint8 word ) 
{
    ASSERT( word > m_lastWord || ( m_dataSize == 0 && m_suffixes.size() == 0 ) );
    Uint8 index = super_type::GetConfig().Word2Index( word );
    if( index > m_lastIndex ) Purge();
    m_suffixes.push_back( super_type::GetConfig().Word2Suffix( word ) );
    m_lastIndex = index;
    m_lastWord = word;
    m_addCount++;
}

template<class TPackConfig>
inline void CPackedBitmaskDataAppender<TPackConfig>::Purge()
{
    ASSERT( m_suffixes.size() <= (super_type::GetConfig().GetSuffixMask() + 1) );
    ASSERT( m_dataSize <= super_type::GetConfig().GetOffsetMask() );
    Uint8 longestList = super_type::GetConfig().GetLongestList();
    Uint8 listSize = m_suffixes.size();
    if( listSize > longestList ) {
        listSize = super_type::GetConfig().GetCountMask();
    }

    CBitHacks::SetBits( super_type::m_table, m_lastIndex*super_type::GetConfig().GetRecordBits(),  super_type::GetConfig().GetCountBits(), listSize );
    CBitHacks::SetBits( super_type::m_table, m_lastIndex*super_type::GetConfig().GetRecordBits() + super_type::GetConfig().GetCountBits(), super_type::GetConfig().GetOffsetBits(), m_dataSize );
    if( m_suffixes.size() == 0 ) return;

    Uint8 extra = super_type::GetConfig().GetDataEntryBytes( m_suffixes.size() );
    super_type::m_data = (Uint1*)realloc( super_type::m_data, m_dataSize + extra );
    memset( super_type::m_data + m_dataSize, 0, extra );
    if( m_suffixes.size() > longestList ) {
        ITERATE( TSuffixes, sfx, m_suffixes ) {
            CBitHacks::SetBits( super_type::m_data + m_dataSize, *sfx, 1, 1 );
        }
    } else {
        for( unsigned i = 0; i < m_suffixes.size(); ++i ) 
            CBitHacks::SetBits( super_type::m_data + m_dataSize, super_type::GetConfig().GetSuffixBits()*i, super_type::GetConfig().GetSuffixBits(), m_suffixes[i] );
    }
    if( m_useHashOffsets && m_suffixes.size() ) {
        Uint8 hash = Uint8( m_suffixes.size() ) << (2*super_type::GetConfig().GetSuffixBits()); 
        hash |= m_suffixes.front() << super_type::GetConfig().GetSuffixBits();
        hash |= m_suffixes.back();
        TOffsetHash::iterator x = m_offsetHash.find( hash );
        if( x == m_offsetHash.end() ) {
            x = m_offsetHash.insert( x, make_pair( hash, TOffsets(1) ) );
            x->second[0] = m_dataSize;
        } else {
            // since length is part of key, they all are of the same length
            // and data block is of the same size
            ITERATE( TOffsets, o, x->second ) {
                if( memcmp( super_type::m_data + *o, super_type::m_data + m_dataSize, extra ) == 0 ) {
                    CBitHacks::SetBits( super_type::m_table, m_lastIndex*super_type::GetConfig().GetRecordBits() + super_type::GetConfig().GetCountBits(), super_type::GetConfig().GetOffsetBits(), *o );
                    m_suffixes.clear();
                    m_savedBytes += extra;
                    return;
                }
            }
            x->second.push_back( m_dataSize );
        }
    }
    m_dataSize += extra;
    m_bitCount += m_suffixes.size();
    m_suffixes.clear();
}

//////////////////////////////////////////////////////////////////////
// CBitmaskFileHeader

inline CBitmaskFileHeader::CBitmaskFileHeader( const CVarBmPackConfig& cfg ) :
    CVarBmPackConfig( cfg ), 
    m_byteOrder( 0x0123456789abcdefULL ),
    m_headerSize(0),
    m_bitsPerUnit(0),
    m_unitsTotal(0),
    m_maxAmb(0),
    m_winStep(0),
    m_wordsPresent(0),
    m_basePattern(0), 
    m_ncbi2naPattern(0), 
    m_winLength(0),
    m_flags(0),
    m_contentOffset(0),
    m_dataSize(0),
    m_signature(36),
    m_signatureLength(0)
{ x_InitSignature(); }

inline void CBitmaskFileHeader::x_InitSignature()
{
    static char sig[] = "oligoFAR.hash.word.bit.mask:0.0.0\0\0\0";
    m_signature.resize( sizeof(sig)-1 );
    copy( sig, sig + sizeof(sig)-1, &m_signature[0] );
    m_signatureLength = strchr( sig, ':' ) - sig;
}

inline void CBitmaskFileHeader::SetFileMinor( int ver )
{
    //DEBUG( DISPLAY( ver ) << DISPLAY( m_winLength ) << DISPLAY( GetWordBases() ) << DISPLAY( GetPrefixBits() ) );
    if( ver < 1 && m_winLength != GetWordBases() ) ver = 1;
    if( ver < 2 && GetPrefixBits() > 0 ) ver = 2;
    ASSERT( ver < 3 );
    m_signature[m_signatureLength + 3] = ver + '0';
    m_byteOrder = ver < 2 ? 0x01234567 : 0x0123456789abcdefULL;
    //DEBUG( DISPLAY( ver ) << DISPLAY( &m_signature[0] ) << DISPLAY( GetFileMinor() ) );
}

inline void CBitmaskFileHeader::SetPattern( Uint8 pattern ) 
{
    if( pattern ) while( (pattern & 1) == 0 ) pattern>>=1;
    m_basePattern = pattern;
    m_winLength = CBitHacks::LargestBitPos( pattern ) + 1;
    int wordSize = CBitHacks::BitCount4( pattern );
    int wordBits = 2*wordSize;
    m_ncbi2naPattern = CBitHacks::DuplicateBits<Uint8>( pattern );
    m_pfxBits = min( m_pfxBits.Get(), Uint2(wordBits) );
    m_sfxBits = wordBits - m_pfxBits.Get();
    //DEBUG( DISPLAYh( m_basePattern ) << DISPLAYh( m_ncbi2naPattern ) ); 
    //m_wordMask = GetWordMask();
}

inline void CBitmaskFileHeader::SetPrefixBits( int bits ) 
{
    int wordBits = GetWordBits();
    m_pfxBits = min( bits, wordBits );
    m_sfxBits = wordBits - m_pfxBits.Get();
}

template <int szof, class T>
inline void CBitmaskFileHeader::x_Read( istream& in, T* dest, int cnt )
{
    ASSERT( sizeof(T) == szof );
    in.read( (char*)dest, sizeof(T)*cnt );
}

template <int szof, class T>
inline void CBitmaskFileHeader::x_Write( ostream& out, const T* src, int cnt )
{
    ASSERT( sizeof(T) == szof );
    out.write( (const char*)src, sizeof(T)*cnt );
}

inline void CBitmaskFileHeader::Read( ifstream& in ) 
{
    in.seekg(0);
    x_InitSignature();
    x_Read<8>( in, &m_byteOrder );
    if( m_byteOrder != 0x0123456789abcdefULL && m_byteOrder != 0x01234567 ) 
        THROW( runtime_error, "Wrong bitmask file byte order" );
    x_Read<4>( in, &m_headerSize );
    x_Read<4>( in, &m_bitsPerUnit );
    x_Read<8>( in, &m_unitsTotal );
    x_Read<4>( in, &m_maxAmb );
    Uint4 wordSize = 0;
    x_Read<4>( in, &wordSize );
    x_Read<4>( in, &m_winStep );
    x_Read<8>( in, &m_wordsPresent );
    m_basePattern = 0;
    if( in.fail() ) THROW( runtime_error, "Failed to read file header" );
    if( m_headerSize > (Uint8)in.tellg() + GetSignatureSize() - sizeof( m_byteOrder ) )
        x_Read<4>( in, &m_basePattern );
    else
        m_basePattern = CBitHacks::WordFootprint<Uint4>( wordSize );
    
    m_contentOffset = m_headerSize + sizeof( m_byteOrder );

    vector<char> signature( GetSignatureSize() );
    x_Read<1>( in, &signature[0], GetSignatureSize() );

    if( strncmp( &signature[0], &m_signature[0], m_signatureLength ) )
        THROW( runtime_error, "Bad file signature" );

    copy( signature.begin(), signature.end(), &m_signature[0] );

    m_sfxBits = wordSize*2;
    if( in.fail() ) THROW( runtime_error, "Failed to read file header" );
    if( GetFileMajor() == 0 && GetFileMinor() == 2 && GetFileRelease() == 0 ) {
        if( m_headerSize < (Uint8)in.tellg() + 12 )
            THROW( runtime_error, "Flags and prefix should follow signature in file format 0.2.0" );
        x_Read<8>( in, &m_flags );
        x_Read<2>( in, &m_pfxBits );
        x_Read<2>( in, &m_cntBits );
        x_Read<2>( in, &m_offBits );
        x_Read<8>( in, &m_dataSize );
        if( in.fail() ) THROW( runtime_error, "Failed to read file header" << strerror( errno ) );
    } else {
        m_flags = 0;
        m_pfxBits = 0;
        m_cntBits = 0;
        m_offBits = 0;
        m_dataSize = (GetWordMask()+1)/8;
    }
    SetPattern( m_basePattern );
    SetPrefixBits( m_pfxBits.Get() );
    Reconsile(GetFileMinor());
}   
    
inline void CBitmaskFileHeader::Reconsile( int minor ) 
{
    //m_bitsPerUnit = 32;
    //m_wordsPresent = m_wordMask + 1;
    if( m_bitsPerUnit == 0 ) THROW( logic_error, "CBitmaskFileHeader::Reconsile(): need bitls per unit to be set" );
    m_winStep = 1;
    m_unitsTotal = (GetWordMask()+1)/m_bitsPerUnit;
    //DEBUG( DISPLAY( minor ) );
    if( minor > 2 ) minor = 2;
    if( minor < 0 ) minor = 0;
    //DEBUG( DISPLAY( minor ) );
    SetFileMinor(minor);
    if( minor < 2 && IsCompressed() )
        THROW( runtime_error, "Can't write compressed bitmask data in compatibility mode" );
    if( IsCompressed() ) {
        m_bitsPerUnit = GetRecordBits();
    } 
    
    m_headerSize = 
        sizeof( m_headerSize ) +
        sizeof( m_bitsPerUnit ) +
        sizeof( m_unitsTotal ) +
        sizeof( m_maxAmb ) +
        sizeof( Uint4 ) + // wordSize
        sizeof( m_winStep ) +
        sizeof( m_wordsPresent );
    //DEBUG( DISPLAY( GetFileMinor() ) );
    ASSERT( GetFileMajor() == 0 && GetFileRelease() == 0 );
    switch( GetFileMinor() ) {
        case 2: m_headerSize += sizeof( m_flags ) + sizeof( Uint2 /*m_prefixBits*/ ) + 
                sizeof( Uint2 /*m_countBits*/ ) + sizeof( Uint2 /*m_tableOffsetBits*/ ) + 
                    sizeof( m_dataSize );
        case 1: m_headerSize += sizeof( m_basePattern );
        case 0: m_headerSize += GetSignatureSize(); 
            break; 
        default: THROW( logic_error, "Unrecognzied requested save file version" );
    }
    m_contentOffset = (m_headerSize + sizeof( m_byteOrder ) + 15)&(~Uint8(0x0f));
}

inline void CBitmaskFileHeader::Write( ofstream& out ) 
{
    out.seekp(0);
    x_Write<8>( out, &m_byteOrder );
    x_Write<4>( out, &m_headerSize );
    x_Write<4>( out, &m_bitsPerUnit );
    x_Write<8>( out, &m_unitsTotal );
    x_Write<4>( out, &m_maxAmb );
    Uint4 wordSize = GetWordBases();
    x_Write<4>( out, &wordSize );
    x_Write<4>( out, &m_winStep );
    x_Write<8>( out, &m_wordsPresent );
    if( GetFileMinor() > 0 ) x_Write<4>( out, &m_basePattern );
    x_Write<1>( out, GetSignature(), GetSignatureSize() );
    if( GetFileMinor() > 1 ) { 
        x_Write<8>( out, &m_flags );
        x_Write<2>( out, &m_pfxBits );
        x_Write<2>( out, &m_cntBits );
        x_Write<2>( out, &m_offBits );
        x_Write<8>( out, &m_dataSize );
    }
    if( out.fail() )
        THROW( runtime_error, "Failed to write bitmask file info" );
    ASSERT( GetContentOffset() >= GetHeaderSize() + sizeof( m_byteOrder ) );
    out.seekp( GetContentOffset() );
}

inline void CBitmaskFileHeader::Print( ostream& out ) 
{
#define PRINT( a ) out << setw( 40 ) << left << #a << ":\t" << a << "\n"
#define PRINTH( a ) out << setw( 40 ) << left << #a << ":\t0x" << hex << a << dec << "\n"
    PRINT( GetPrefixBits() );
    PRINT( GetSuffixBits() );
    PRINT( GetOffsetBits() );
    PRINT( GetCountBits() );
    PRINT( GetWordBits() );
    PRINT( GetRecordBits() );
    PRINT( GetRecordCount() );
    PRINT( GetTableBytesPacked() );
    PRINT( GetTableBytesAligned() );
    PRINT( GetLongestList() );
    PRINT( GetBitmaskBits() );
    PRINT( GetBitmaskBytes() );
    PRINT( GetWordBases() );
    PRINTH( GetWordBaseMask() );
    PRINTH( GetWordMask() );
    PRINTH( GetIndexMask() );
    PRINTH( GetSuffixMask() );
    PRINTH( GetPrefixMask() );
    PRINT( GetSignatureLength() );
    PRINT( GetSignatureSize() );
    PRINT( GetSignature() );
    PRINT( GetFileMajor() );
    PRINT( GetFileMinor() );
    PRINT( GetFileRelease() );
    PRINTH( GetByteOrder() );
    PRINT( GetHeaderSize() );
    PRINT( GetBitsPerUnit() );
    PRINT( GetUnitsTotal() );
    PRINT( GetMaxAmb() );
    PRINT( GetWinStep() );
    PRINT( GetWordsPresent() );
    PRINTH( GetBasePattern() );
    PRINTH( GetNcbi2naPattern() );
    PRINT( GetWindowLength() );
    PRINT( GetFlags() );
    PRINT( GetContentOffset() );
    PRINT( GetTableOffset() );
    PRINT( GetDataOffset() );
    PRINT( GetDataSize() );
    if( IsCompressed() ) {
        PRINT( GetDataEntryBytes( 1 ) );
        PRINT( GetDataEntryBytes( 2 ) );
        PRINT( GetDataEntryBytes( 4 ) );
        PRINT( GetDataEntryBytes( GetLongestList() ) );
        PRINT( GetDataEntryBytes( GetLongestList() + 1 ) );
        PRINT( GetDataEntryBytes( GetSuffixMask() ) );
        PRINT( GetDataEntryBytes( GetSuffixMask() + 1 ) );
    }
#undef PRINT
#undef PRINTH
}

inline CBitmaskAccess::CBitmaskAccess() : 
    m_flatBm(0,0), m_packedBm( CVarBmPackConfig( 0, 0, 0, 0 ), 0, 0 ),
    m_xflags( 0 )
{}

inline void CBitmaskAccess::ReleaseFile() 
{
    if( IsCompressed() ) {
        delete[] (const char*)m_packedBm.GetTable();
        delete[] (const char*)m_packedBm.GetData();
        m_packedBm = TPackedBitmap( CVarBmPackConfig( 0,0,0,0 ), 0,0 );
    } else {
        if( m_xflags & fUseMmap ) {
            munmap( ((char*)m_flatBm.GetData()) - GetDataOffset(), GetDataSize() + GetDataOffset() );
            m_xflags &= ~fUseMmap;
        } else {
            delete[] m_flatBm.GetData();
        }
        m_flatBm = TFlatBitmap( 0, 0 );
    }
}

inline void CBitmaskAccess::AttachFile( const string& fname, int flags, IProgressIndicator * p )
{
    ReleaseFile();
    ifstream in( fname.c_str() );
    Read( in );
    //Print( cerr << "AFTER Read()\n" );
    m_flags = flags;
    if( IsCompressed() ) {
        if( flags & fUseMmap ) flags &= ~fUseMmap;
        if( flags & fUseMmap ) {
            THROW( logic_error, "Mmap of compressed files is not implemented yet" );
        } else {
            char * table = new char[ GetTableBytesAligned() ];
            char * data = new char[ GetDataSize() ];
            in.seekg( GetTableOffset() );
            in.read( table, GetTableBytesAligned() );
            if( (Uint8)in.gcount() != GetTableBytesAligned() )
                THROW( runtime_error, "Failed to read " << GetTableBytesAligned() << " bytes bitmask index from file " << fname << ": " << strerror( errno ) );
            in.seekg( GetDataOffset() );
#ifndef USE_LARGE_READ
#define USE_LARGE_READ 0
#endif
#if USE_LARGE_READ
            in.read( data, GetDataSize() );
            if( (Uint8)in.gcount() != GetDataSize() )
                THROW( runtime_error, "Failed to read " << GetDataSize() << " bytes bitmask data from file " << fname << ": " << strerror( errno ) );
#else
            const size_t bsize = 65536;
            size_t off = 0;
            while( off < GetDataSize() ) {
                size_t sz = min( bsize, GetDataSize() - off );
                in.read( data + off, sz );
                if( (size_t)in.gcount() != sz ) 
                    THROW( runtime_error, "Failed to read " << sz << " bytes bitmask data from data block offset " << off << " of file " << fname << ": " << strerror( errno ) );
                off += sz;
            }
#endif

            m_packedBm = TPackedBitmap( *this, table, data );
            m_xflags = flags;
        }
    } else {
        if( flags & fUseMmap ) { 
           int fd = open( fname.c_str(), O_RDONLY );
           if( fd == -1 ) 
               THROW( runtime_error, "Oops... failed to open() file " << fname << ": " << strerror( errno ) );
           Uint8 len = TFlatBitmap::GetBytes( GetWordBases() ) + GetDataOffset();
           char * d = (char*)mmap( 0, len, PROT_READ, MAP_PRIVATE|MAP_NORESERVE, fd, 0 );
           if( d == MAP_FAILED ) {
               close( fd );
               THROW( runtime_error, "Oops... failed to mmap() file " << fname << ": " << strerror( errno ) );
           }
           madvise( d, len, (flags&fMmapSequential)?(MADV_SEQUENTIAL|MADV_DONTNEED):(MADV_RANDOM|MADV_WILLNEED) );
           close( fd );
           m_xflags|=fUseMmap;
           m_flatBm = TFlatBitmap( d + GetDataOffset(), GetWordBases() );
           //DEBUG( "Data (" << TFlatBitmap::GetBytes( GetWordBases() ) << " bytes) is memory mapped to " << (void*)d );
        } else {
           TFlatBitmap::data_type * data = TFlatBitmap::Allocate( GetWordBases() );
           //DEBUG( "Data (" << TFlatBitmap::GetBytes( GetWordBases() ) << " bytes) allocated at " << data );
           m_flatBm = TFlatBitmap( data, GetWordBases() );
           in.seekg( GetDataOffset() );
           for( Uint8 off = 0; off < m_flatBm.GetDataBytes(); ) {
               in.read( ((char*)data) + off, min( m_flatBm.GetDataBytes() - off, Uint8(1024*1024) ) );
               if( in.gcount() == 0 ) THROW( runtime_error, "Failed to read data from " << fname << ": " << strerror( errno ) );
               off += in.gcount();
               if( p ) p->Increment();
           }
        }

        /*
        for( Uint8 b = 0; b < bytes/sizeof( TFlatBitmap::data_type); ) {
            for( int k = 0; k < 8; ++k, ++b ) {
                cerr << " " << hex << setw(sizeof( TFlatBitmap::data_type )*2) << setfill('0') << data[b] << dec;
            }
            cerr << "\n";
        }
        */
    }
}

template<class TData, class TCfg>
inline void CBitmaskWriter::InitHeader( CBitmaskFileHeader& header, const TData& data, const TCfg& cfg )
{
    CCounter counter;
    data.ForEachSet( counter );
    
    if( m_basePattern == 0 ) m_basePattern = CBitHacks::WordFootprint<Uint8>( cfg.GetWordBaseMask() );
    m_ncbi2naPattern = CBitHacks::DuplicateBits<Uint8>( m_basePattern );
    //DEBUG( DISPLAY( cfg.GetWordBases() ) << hex << DISPLAY( m_basePattern ) << DISPLAY( cfg.GetWordBaseMask() ) << dec );
    if( CBitHacks::BitCount8( m_basePattern ) > (int)cfg.GetWordBases() ) header.SetFileMinor(1);
    else if( CBitHacks::BitCount8( m_basePattern ) < (int)cfg.GetWordBases() ) THROW( logic_error, "Pattern is smaller then word!" );
    else header.SetFileMinor(0);
    header.SetPattern( m_basePattern );
    header.SetWordsPresent( counter.GetCount() );
}

template<class Unit>
inline void CBitmaskWriter::Write( const string& fname, const CFlatBitmaskData<Unit>& data, IProgressIndicator * p )
{
    CBitmaskFileHeader header;
    InitHeader( header, data );

    header.SetPrefixBits(0);
    header.SetUnitBits(sizeof(Unit)*8);
    //DEBUG( "Reconsiling to " << DISPLAY( m_minor ) );
    header.Reconsile(m_minor);

    ofstream out( fname.c_str() );
    //header.Print( cerr << "To Write:\n" );
    header.Write( out );
    out.seekp( header.GetDataOffset() );
    for( Uint8 off = 0; off < data.GetDataBytes();  ) {
        Uint8 sz = min( Uint8(1024*1024), data.GetDataBytes() - off );
        out.write(((const char*) data.GetData()) + off, sz );
        if( out.fail() ) THROW( runtime_error, "Failed to write data to file " << fname << ": " << strerror( errno ) );
        off += sz;
        if( p ) p->Increment();
    }
}

template<class Config>
inline void CBitmaskWriter::Write( const string& fname, const CPackedBitmaskData<Config>& data, Uint8 dsize, IProgressIndicator * p )
{
    CBitmaskFileHeader header( data.GetConfig() );
    CPackedBitmaskDataReader<Config> reader( data.GetConfig(), data.GetTable(), data.GetData() );
    InitHeader( header, reader, data.GetConfig() );

    header.SetUnitBits(8);
    header.SetPrefixBits(data.GetConfig().GetPrefixBits());
    header.Reconsile(2);
    header.SetDataSize( dsize );

    ofstream out( fname.c_str() );
    //header.Print( cerr << "To Write:\n" );
    header.Write( out );
    out.seekp( header.GetTableOffset() );
    out.write((const char*) data.GetTable(), data.GetConfig().GetTableBytesPacked() );
    out.seekp( header.GetDataOffset() );
    out.write((const char*) data.GetData(), dsize );
}

inline bool CBitmaskAccessExt::HasWord( Uint8 word ) const
{
    if( m_rqWordSize == 0 || m_rqWordSize == GetWordBases() ) {
        for( Uint4 i = GetWordBases(); i <= m_rqWordSize; ++i, (word >>= 2) )
            if( !HasExactWord( word ) ) { return false; }
        return true; // all subwords exist
    } else {
        int k = GetWordBases() - m_rqWordSize;
        Uint8 wstart = word << (2*k);
        Uint8 wend = (word + 1) << (2*k);
        for( Uint8 w = wstart; w != wend; ++w ) 
            if( HasExactWord( w ) ) { return true; }
        return false;
    }
}

END_OLIGOFAR_SCOPES

#endif
