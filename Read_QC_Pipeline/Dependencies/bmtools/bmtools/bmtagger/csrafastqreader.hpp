#ifndef OLIGOFAR_CSRAFASTQREADER__HPP
#define OLIGOFAR_CSRAFASTQREADER__HPP

/*  $Id: csrafastqreader.hpp 382804 2012-12-10 02:38:55Z rotmistr $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Kirill Rotmistrovsky
 *
 *
 */

#ifndef WITH_SRA
#define WITH_SRA 0
#endif
#if WITH_SRA

#ifndef SRA_VERSION
#define SRA_VERSION 2
#endif


#include <klib/defs.h> 
#include <vdb/types.h> 
#include <klib/rc.h> 
#if SRA_VERSION == 1
#include <sra/ncbi-sradb.h>
#else
#define SRAErrMake( a ) (a)
#define SRAMgrOpenRunRead( a, b, c ) SRAMgrOpenTableRead( (a), (b), (c) )
#include <klib/writer.h>
#endif
#include <sra/sradb.h>
#include <sra/fastq.h>

//#include <netinet/in.h>
//#include <util/range.hpp>
//#include <util/value_convert.hpp>

#include <vector>
#include <list>
#include <set>

#include "of-debug.hpp"
#include "of-scopes.hpp"

BEGIN_OLIGOFAR_SCOPES

class CSRAFastqReader
{
public:
    // temporary - sradb does not have readfilter_redacted yet
    typedef pair<short,short> TRange;
    enum EReadFilter { eReadFilter_pass, eReadFilter_reject, eReadFilter_criteria, eReadFilter_redacted }; 
    enum EFlags { 
        fReadId_includeAccession = 0x01, 
        fReadDb_useFastqReader = 0x02, 
        fReadDb_readClipInfo = 0x04,
        fReadDb_needQualities = 0x08,
        fFlags_default = fReadId_includeAccession 
    };

    CSRAFastqReader( const string& accession, Uint4 flags = fFlags_default, int shortestRead = 25 );
    ~CSRAFastqReader();
    rc_t GetRc() const { return m_rc; }
    int GetError() const { return SRAErrMake( GetRc() ); }
#if SRA_VERSION == 1
    const char * GetErrorMsg() const { return SRAErrToEnglish( GetError(), 0 ); }
#else
    string GetErrorMsg() const { 
        char buff[4096];
        size_t n = 0;
        RCExplain( GetError(), buff, sizeof(buff), &n );
        return string( buff, n );
    }
#endif
    spotid_t GetMaxSpotId() const { return m_maxSpotId; }
    spotid_t GetSpotId() const { return m_spotId; }
    const string& GetReadId() const { return m_readId; }
    const string& GetRead( int c ) const { return m_read[c]; }
    const string& GetQual( int c ) const { return m_qual[c]; }
    const TRange& GetClip( int c ) const { ASSERT( c >= 0 && c < (int)m_read.size() ); return m_clip[c]; }
    int GetComponents() const { return m_read.size(); }
    bool Exists() const { return m_exists; }
    bool Loaded() const { return m_loaded; }
    bool NextRead();
    bool FetchRead();
    void Rewind();
    void SetRqQuality( bool c ) { m_rqQuality = c; }
    template<class Iterator>
    void SetSpotIdList( Iterator b, Iterator e ) {
        m_spotIdList.clear();
        while( b != e ) m_spotIdList.insert( *b++ ); // Convert<Uint8>( *b++ ) );
    }

    void SetFlags( Uint4 flags, bool on ) { if( on ) m_flags |= flags; else m_flags &= ~flags; }
    void AssignFlags( Uint4 flags ) { m_flags = flags; }
    Uint4 GetFlags() const { return m_flags; }

    EReadFilter GetReadFilterValue() const { return m_rfValue; }
protected:
    string m_accession;
    SRAMgr const * m_sraMgr;
    SRATable const * m_sraTable;
    SRAColumn const * m_readSeq;
    SRAColumn const * m_readLen;
    SRAColumn const * m_readType;
    SRAColumn const * m_readFilter;
#if SRA_VERSION == 1 
    SRAColumn const * m_clipLeft;
    SRAColumn const * m_clipRight;
#else
    SRAColumn const * m_trimStart;
    SRAColumn const * m_trimLen;
#endif
    FastqReader const * m_fastqReader;
    rc_t m_rc;
    spotid_t m_maxSpotId;
    spotid_t m_spotId;
    vector<string> m_read;
    vector<string> m_qual;
    vector<TRange> m_clip;
    string m_readId;
    Uint4 m_flags;
    bool m_rqQuality;
    bool m_exists;
    bool m_loaded;
    EReadFilter m_rfValue;
    set<Uint8> m_spotIdList;
};

inline CSRAFastqReader::CSRAFastqReader( const string& accession, Uint4 flags, int shortestRead ) :
    m_accession( accession ),
    m_sraMgr( 0 ),
    m_sraTable( 0 ),
    m_readSeq( 0 ),
    m_readLen( 0 ),
    m_readType( 0 ),
    m_readFilter( 0 ),
#if SRA_VERSION == 1
    m_clipLeft( 0 ),
    m_clipRight( 0 ),
#else
    m_trimStart( 0 ),
    m_trimLen( 0 ),
#endif
    m_fastqReader( 0 ),
    m_rc( 0 ),
    m_maxSpotId( 0 ),
    m_spotId( 0 ),
    m_flags( flags ),
    m_rqQuality( false ),
    m_exists( false ),
    m_loaded( false )
{
    m_rc = SRAMgrMakeRead( &m_sraMgr );
    if( m_rc ) THROW( runtime_error, "Failed to init SRA read: " << GetErrorMsg() );
    m_rc = SRAMgrOpenRunRead( m_sraMgr, &m_sraTable, m_accession.c_str() );
    if( m_rc ) THROW( runtime_error, "Failed to open SRA run " << m_accession << ": " << GetErrorMsg() );
    m_rc = SRATableMaxSpotId( m_sraTable, &m_maxSpotId );
    if( m_rc ) THROW( runtime_error, "Failed to get max spot ID for run " << m_accession << ": " << GetErrorMsg() );
    if( m_flags & fReadDb_useFastqReader ) {
        m_rc = FastqReaderMake( &m_fastqReader, m_sraTable, accession.c_str(), 
                                false, // colorspace
                                false, // original format - keep SRR accession in defline
                                m_flags&fReadDb_needQualities,  // qualities
                                false, // print read label
                                true,  // print read id
                                false, // no_clip = false
                                shortestRead, // even this short
                                33,    // quality conversion offset
                                '\0',  // no color conversion
                                0, 0 );// min and max spot ID - use run info
        if( m_rc ) THROW( runtime_error, "Failed to init FASTQ reader for " << m_accession << ": " << GetErrorMsg() );
        else cerr << "INFO: opened FASTQ reader\n";
    } else {
        // 
        m_rc = SRATableOpenColumnRead( m_sraTable, &m_readSeq, "READ", "INSDC:dna:text" ); 
        if( m_rc ) THROW( runtime_error, "Failed to init READ column of type INSDC:dna:text for " << m_accession << ": " << GetErrorMsg() );
        else cerr << "INFO: opened column READ\n";
        m_rc = SRATableOpenColumnRead( m_sraTable, &m_readLen, "READ_LEN", "INSDC:coord:len" ); // Changed according to schema vdb_uint16_t );
        if( m_rc ) THROW( runtime_error, "Failed to init READ_LEN column for " << m_accession << ": " << GetErrorMsg() );
        else cerr << "INFO: opened column READ_LEN\n";
        m_rc = SRATableOpenColumnRead( m_sraTable, &m_readType, "READ_TYPE", "INSDC:SRA:read_type" ); // sra_read_type_t );
        if( m_rc ) THROW( runtime_error, "Failed to init READ_TYPE column for " << m_accession << ": " << GetErrorMsg() );
        else cerr << "INFO: opened column READ_TYPE\n";
        if( m_flags & fReadDb_readClipInfo ) {
#if SRA_VERSION == 1
            m_rc = SRATableOpenColumnRead( m_sraTable, &m_clipLeft, "CLIP_QUALITY_LEFT", "INSDC:coord:one" );  // vdb_uint16_t );
            if( m_rc ) cerr << "WARNING: Failed to init CLIP_QUALITY_LEFT column for " << m_accession << ": " << GetErrorMsg() << "\n";
            else { cerr << "INFO: opened column CLIP_QUALITY_LEFT\n"; ASSERT( m_clipLeft ); }
            m_rc = SRATableOpenColumnRead( m_sraTable, &m_clipRight, "CLIP_QUALITY_RIGHT", "INSDC:coord:one" ); // vdb_uint16_t );
            if( m_rc ) cerr << "WARNING: Failed to init CLIP_QUALITY_RIGHT column for " << m_accession << ": " << GetErrorMsg() << "\n";
            else { cerr << "INFO: opened column CLIP_QUALITY_RIGHT\n"; ASSERT( m_clipRight ); }
#else 
            m_rc = SRATableOpenColumnRead( m_sraTable, &m_trimStart, "TRIM_START", "INSDC:coord:zero" );  // vdb_uint16_t );
            if( m_rc ) THROW( runtime_error, "Failed to init TRIM_START column for " << m_accession << ": " << GetErrorMsg() );
            else cerr << "INFO: opened column TRIM_START\n";
            m_rc = SRATableOpenColumnRead( m_sraTable, &m_trimLen, "TRIM_LEN", "INSDC:coord:len" ); // vdb_uint16_t );
            if( m_rc ) THROW( runtime_error, "Failed to init TRIM_LEN column for " << m_accession << ": " << GetErrorMsg() );
            else cerr << "INFO: opened column TRIM_LEN\n";
#endif
        }
    }
    m_rc = SRATableOpenColumnRead( m_sraTable, &m_readFilter, "READ_FILTER", "INSDC:SRA:read_filter" ); //sra_read_filter_t );
    if( m_rc ) cerr << "WARNING: Failed to init READ_FILTER column for " << m_accession << ": " << GetErrorMsg() << "\n";
    else cerr << "INFO: opened column READ_FILTER\n";
    m_rc = 0;
}

inline CSRAFastqReader::~CSRAFastqReader() 
{
    FastqReaderWhack( m_fastqReader );
#if SRA_VERSION == 1
    SRAColumnRelease( m_clipRight );
    SRAColumnRelease( m_clipLeft );
#else
    SRAColumnRelease( m_trimLen );
    SRAColumnRelease( m_trimStart );
#endif
    SRAColumnRelease( m_readLen );
    SRAColumnRelease( m_readType );
    SRAColumnRelease( m_readFilter );
    SRATableRelease( m_sraTable );
    SRAMgrRelease( m_sraMgr );
}

inline void CSRAFastqReader::Rewind()
{
    m_loaded = false;
    if( m_fastqReader ) {
        m_rc = FastqReaderFirstSpot( m_fastqReader );
        if( m_rc ) THROW( runtime_error, "Failed to seek to first spot of " << m_accession << ": " << GetErrorMsg() );
    } else {
        if( m_spotIdList.size() ) m_spotId = *m_spotIdList.begin();
        else m_spotId = 1;
    }
    m_exists = true;
}

inline bool CSRAFastqReader::FetchRead()
{
    m_read.clear();
    m_qual.clear();
    m_clip.clear();
    m_readId.clear();

    if( !m_exists ) THROW( runtime_error, "Can't fetch non-existing spot" );

    if( m_fastqReader ) {
        m_rc = FastqReaderCurrentSpot( m_fastqReader, &m_spotId );
        if( m_rc ) THROW( runtime_error, "Failed to get current spot ID for " << m_accession << ": " << GetErrorMsg() );
    } else {}

    m_readId = ( m_flags & fReadId_includeAccession ? m_accession + "." : string("") ) + ToStr( m_spotId );
    
    if( m_readFilter ) {
        INSDC_SRA_read_filter * rfData = 0;
        // uint8_t const * rfData = 0;
        bitsz_t off = 0;
        bitsz_t sz = 0;
        m_rc = SRAColumnRead( m_readFilter, m_spotId, (const void**)&rfData, &off, &sz );
        switch( *rfData ) {
            default: THROW( runtime_error, "Unexpected value " << *rfData << " (off=" << off << ",sz=" << sz << ") for READ_FILTER column of read " << m_spotId << " for accession " << m_accession );
            case eReadFilter_pass:
            case eReadFilter_reject:
            case eReadFilter_criteria:
            case eReadFilter_redacted:
                m_rfValue = EReadFilter( *rfData );
        }
    } else m_rfValue = eReadFilter_pass;

    uint32_t numReads = 0;
    if( m_fastqReader ) {

#if SRA_VERSION==1
        m_rc = FastqReader_SpotInfo( m_fastqReader, NULL, NULL, &numReads );
#else
        uint32_t spotLength = 0;
        m_rc = FastqReader_SpotInfo( m_fastqReader, NULL, NULL, NULL, NULL, &spotLength, &numReads );
#endif
        if( m_rc ) THROW( runtime_error, "Failed to get spot info for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
        typedef list<pair<int,int> > TReadInfo;
        TReadInfo readInfo;
        for( unsigned readId = 1; readId <= numReads; ++readId ) {
            SRAReadTypes readType = SRA_READ_TYPE_TECHNICAL;
#if SRA_VERSION==1
            uint16_t readLength = 0;
            m_rc = FastqReader_SpotReadInfo( m_fastqReader, readId, NULL, &readType, NULL, &readLength, NULL );
#else
            INSDC_coord_len  readLength = 0;
            INSDC_coord_zero readStart = 0;
            m_rc = FastqReader_SpotReadInfo( m_fastqReader, readId, &readType, NULL, NULL, &readStart, &readLength );
#endif
            if( m_rc ) THROW( runtime_error, "Failed to get spot read info for spot " << m_spotId << " read " << readId << " of " << m_accession << ": " << GetErrorMsg() );
            if( readType == SRA_READ_TYPE_BIOLOGICAL ) {
                readInfo.push_back( make_pair( readId, readLength ) );
            }
        }
        ITERATE( TReadInfo, r, readInfo ) {
            m_read.push_back( "" );
            m_qual.push_back( "" );
            m_clip.push_back( TRange() );
            if( r->second == 0 ) continue;
            char read[20000], qual[20000];
            m_rc = FastqReaderBase( m_fastqReader, r->first, read, sizeof( read ), NULL );
            if( m_rc ) THROW( runtime_error, "Failed to get bases for spot " << m_spotId << " read " << r->first << " of " << m_accession << ": " << GetErrorMsg() );
            m_read.back().assign( read, r->second );
            m_clip.back() = TRange( 0, r->second );
            if( m_rqQuality ) {
                FastqReaderQuality( m_fastqReader, r->first, qual, sizeof( qual ), NULL );
                m_qual.back().assign( qual, r->second );
            }
        }
    } else {
        INSDC_dna_text * segBase = 0;
        INSDC_coord_len * segLen = 0;
        INSDC_SRA_read_type * segType = 0;
#if SRA_VERSION == 1
        INSDC_coord_zero * clipLeft  = 0;
        INSDC_coord_zero * clipRight = 0;
#else
        INSDC_coord_zero * trimStart = 0;
        INSDC_coord_len *  trimLen   = 0;
#endif
        bitsz_t off = 0;
        bitsz_t sz = 0;
        m_rc = SRAColumnRead( m_readLen, m_spotId, (const void**)&segLen, &off, &sz );
        if( m_rc ) THROW( runtime_error, "Failed to get segment lengths for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
        ASSERT( off == 0 ); 
        ASSERT( sz % (8*sizeof(INSDC_coord_len)) == 0 );
        numReads = sz/(8*sizeof(INSDC_coord_len));
        m_rc = SRAColumnRead( m_readType, m_spotId, (const void**)&segType, &off, &sz );
        if( m_rc ) THROW( runtime_error, "Failed to get segment types for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
        ASSERT( off == 0 ); 
        ASSERT( sz % (8*sizeof(INSDC_SRA_read_type)) == 0 ); 
        ASSERT( sz/(8*sizeof(INSDC_SRA_read_type)) == numReads );
        m_rc = SRAColumnRead( m_readSeq, m_spotId, (const void**)&segBase, &off, &sz );
        if( m_rc ) THROW( runtime_error, "Failed to get FASTA for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
        ASSERT( off == 0 ); 
        ASSERT( sz % (8*sizeof(INSDC_dna_text)) == 0 ); 
        int wholeLen = sz/8;
        int cleft = 0;
        int cright = wholeLen;
#if SRA_VERSION == 1
        if( m_clipLeft && m_clipRight ) {
            m_rc = SRAColumnRead( m_clipRight, m_spotId, (const void**)&clipRight, &off, &sz );
            if( m_rc ) THROW( runtime_error, "Failed to get CLIP_QUALITY_RIGHT for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
            ASSERT( off == 0 ); ASSERT( sz % 16 == 0 );
            cright = *clipRight; // 1-based, cright is the first bad base
            m_rc = SRAColumnRead( m_clipLeft, m_spotId, (const void**)&clipLeft, &off, &sz );
            if( m_rc ) THROW( runtime_error, "Failed to get CLIP_QUALITY_LEFT for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
            ASSERT( off == 0 ); ASSERT( sz % 16 == 0 );
            cleft = *clipLeft - 1; // 1-based, cleft is the first good base
        }
#else
        if( m_flags & fReadDb_readClipInfo ) {
            m_rc = SRAColumnRead( m_trimStart, m_spotId, (const void**)&trimStart, &off, &sz );
            if( m_rc ) THROW( runtime_error, "Failed to get TRIM_START for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
            ASSERT( off == 0 ); ASSERT( sz % 16 == 0 );
            m_rc = SRAColumnRead( m_trimLen, m_spotId, (const void**)&trimLen, &off, &sz );
            if( m_rc ) THROW( runtime_error, "Failed to get TRIM_LEN for spot " << m_spotId << " of " << m_accession << ": " << GetErrorMsg() );
            ASSERT( off == 0 ); ASSERT( sz % 16 == 0 );
            cleft = *trimStart; 
            cright = *trimStart + *trimLen;
            ASSERT( cleft >= 0 );
            ASSERT( cright >= cleft );
            ASSERT( cright <= wholeLen );
        }
#endif
        for( unsigned i = 0, off = 0; i < numReads; ++i ) {
            // DEBUG( DISPLAY( m_spotId ) << DISPLAY( wholeLen ) << DISPLAY( i ) << DISPLAY( segLen[i] ) << DISPLAY( int(segType[i]) ) << DISPLAY( off ) << DISPLAY( wholeLen - off ) );
            int slen = segLen[i];
            ASSERT( slen + (int)off <= wholeLen );
            //const int type_mask = SRA_READ_TYPE_BIOLOGICAL|SRA_READ_TYPE_TECHNICAL;
            const int strand_mask = SRA_READ_TYPE_FORWARD|SRA_READ_TYPE_REVERSE;
            switch( segType[i]&(~strand_mask) ) {
            case SRA_READ_TYPE_BIOLOGICAL: 
                m_read.push_back(""); m_qual.push_back(""); m_clip.push_back( TRange( 0, slen ) );
                //if( clipRight && clipLeft ) {
                if( cright < (int)off ) m_clip.back().second = 0;
                else m_clip.back().second = min( slen, cright - (int)off );
                if( cleft > (int)off + slen ) m_clip.back().first = slen;
                else m_clip.back().first = max( 0, cleft - (int)off );
                //}
#if 0           // SRA stores sequence is the native sequencer order which is an expectation for fastq
                if((segType[i]&strand_mask) == SRA_READ_TYPE_REVERSE) {
                    static char rctable[] = 
                        "................""................"
                        ".............-..""................"
                        ".TVGH..CD..M.KN.""..YSAABW.R......"
                        ".tvgh..cd..m.kn.""..ysaabw.r......"
                        "................""................"
                        "................""................"
                        "................""................"
                        "................""................";
                    m_read.back().reserve( slen );
                    for( size_t x = off, X = off+slen; x < X; ++x ) 
                        m_read.back().push_back( rctable[(unsigned char)segBase[x]] );
                } else 
#endif
                    m_read.back().assign( segBase + off, slen );  
                //cerr << "DEBUG:" << DISPLAY( m_read.size() ) << DISPLAY( m_read.back() ) 
                //     << DISPLAY( m_clip.back().first ) << DISPLAY( m_clip.back().second ) << "\n";
                break;
            case SRA_READ_TYPE_TECHNICAL: break;
            default: THROW( logic_error, "Unexpected segment type value of " << (int)segType[i] << " for segment " << i << " spot " << m_spotId << " of " << m_accession );
            }
            off += slen;
        }
    }

    return m_loaded = (m_rc == 0);
}

inline bool CSRAFastqReader::NextRead()
{
    m_read.clear();
    m_qual.clear();
    m_clip.clear();
    m_readId.clear();
    
    m_loaded = false;
    if( m_fastqReader ) {
        m_rc = FastqReaderNextSpot( m_fastqReader );
        if( GetRCObject( m_rc ) == rcRow && GetRCState( m_rc ) == rcExhausted ) { m_rc = 0; return m_exists = false; }
        if( m_rc ) THROW( runtime_error, "Failed to fetch next spot of " << m_accession << ": " << GetErrorMsg() );
    } else {
        if( m_spotIdList.size() ) {
            set<Uint8>::const_iterator x = m_spotIdList.find( m_spotId );
            if( x == m_spotIdList.end() ) return m_exists = false;
            x++;
            if( x == m_spotIdList.end() ) return m_exists = false;
            m_spotId = *x;
        } else m_spotId++;
        if( m_spotId > m_maxSpotId ) return m_exists = false;
    }
    return m_exists = true;
}

END_OLIGOFAR_SCOPES

#endif

#endif

