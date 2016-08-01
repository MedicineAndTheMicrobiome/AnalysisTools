#include "cshortreader.hpp"
#include <fstream>

USING_OLIGOFAR_SCOPES;

////////////////////////////////////////////////////////////////////////
// CPackedFile

CPackedFile::CPackedFile( const string& name ) : m_name( name )
{}

istream& CPackedFile::SetStream() 
{
    if( m_plain.get() == 0 ) {
        m_plain.reset( new ifstream( m_name.c_str() ) );
#ifndef STANDALONE
        if( m_name.length() > 4 && m_name.substr( m_name.length() - 4 ) == ".bz2" )
            m_packed.reset( new CDecompressIStream( *m_plain, CCompressStream::eBZip2 ) );
        else if( m_name.length() > 3 && m_name.substr( m_name.length() - 3 ) == ".gz" )
            m_packed.reset( new CDecompressIStream( *m_plain, CCompressStream::eGZipFile ) );
#endif
    }
    return GetStream();
}

istream& CPackedFile::GetStream() const
{
#ifndef STANDALONE
    if( m_packed.get() ) return *m_packed;
#endif
    if( m_plain.get() ) return *m_plain;
    THROW( logic_error, "Failed to create input stream for " + m_name );
}

////////////////////////////////////////////////////////////////////////
// IShortReader

void IShortReader::SetStream( TFileStream& stream, const string& name )
{
    stream.reset( new CPackedFile( name ) ); // ifstream( name.c_str() ) );
    if( stream->SetStream().fail() )
        THROW( runtime_error, "Error: failed to open file " << name << " for reading" );
}

////////////////////////////////////////////////////////////////////////
// CReaderFactory

unsigned CReaderFactory::x_MkFlags() const
{
    unsigned flags = 0;
    switch( m_qualityChannels ) {
    case 0: flags |= fQuality_0; break;
    case 1: flags |= fQuality_1; break;
    case 4: flags |= fQuality_4; break;
    default: THROW( runtime_error, "Invalid number of quality channels specified: " << m_qualityChannels );
    }
    if( m_colorSpace ) flags |= fColorspace;
    if( m_readIdFile.size() ) flags |= fReadIdFileSet;
    if( m_readDataFile1.size() ) flags |= fReadDataFile1;
    if( m_readDataFile2.size() ) flags |= fReadDataFile2;
    return flags;
}

IShortReader * CReaderFactory::CreateReader()
{
    switch( x_MkFlags() ) {
        case fQuality_0|fReadIdFileSet:
        case fQuality_1|fReadIdFileSet:
        case fQuality_0|fReadIdFileSet|fColorspace:
            cerr << "* Notice: creating CColFileReader( " << m_qualityChannels << ", " << m_colorSpace << ", " << m_readIdFile << ")\n";
            return new CColFileReader( (bool)m_qualityChannels, m_colorSpace, m_readIdFile );
        case fQuality_0|fReadDataFile1:
        case fQuality_1|fReadDataFile1:
        case fQuality_0|fReadDataFile1|fReadDataFile2:
        case fQuality_1|fReadDataFile1|fReadDataFile2:
        case fQuality_0|fReadDataFile1|fColorspace:
        case fQuality_0|fReadDataFile1|fReadDataFile2|fColorspace:
            cerr << "* Notice: creating CFastqFileReader( " << m_qualityChannels << ", " << m_colorSpace << ", " << m_readDataFile1 << ", " << m_readDataFile2 << ")\n";
            return new CFastqFileReader( (bool)m_qualityChannels, m_colorSpace, m_readDataFile1, m_readDataFile2 );
        case fQuality_4|fReadIdFileSet|fReadDataFile1:
        case fQuality_4|fReadIdFileSet|fReadDataFile1|fReadDataFile2:
            cerr << "* Notice: creating CSolexaFileReader( " << m_readIdFile << ", " << m_readDataFile1 << ", " << m_readDataFile2 << ")\n";
            return new CSolexaFileReader( m_readIdFile, m_readDataFile1, m_readDataFile2 );
        default:
            return 0;
    }
}

////////////////////////////////////////////////////////////////////////
// CColFileReader

CColFileReader::CColFileReader( bool quality, bool colorspace, const string& filename )
{
    if( colorspace ) {
        if( quality ) THROW( logic_error, "Oops... don't know how to use quality scores in colorspace" );
        m_coding = CSeqCoding::eCoding_colorsp;
    } else m_coding = quality ? CSeqCoding::eCoding_ncbiqna : CSeqCoding::eCoding_ncbi8na;
    
    SetStream( m_input, filename );
}

bool CColFileReader::NextRead()
{
    m_readId = m_readData[0] = m_readData[1] = "";
    string buff;
    while( getline( m_input->SetStream(), buff ) ) {
        buff = Trim( buff );
        if( buff.length() == 0 || buff[0] == '#' ) continue;
        istringstream in( buff );
        in >> m_readId >> m_readData[0];
        if( !in.eof() ) in >> m_readData[1];
        if( m_readData[1] == "-" ) m_readData[1] = "";
        if( m_coding == CSeqCoding::eCoding_ncbiqna ) {
            string qf, qr;
            in >> qf;
            if( !in.eof() ) in >> qr;
            if( qr == "-" ) qr = "";
            if( m_readData[0].length() != qf.length() )
                THROW( runtime_error, "Error: quality data for read " << m_readId << " (1) does not correspond to sequence length" );
            if( m_readData[1].length() != qr.length() )
                THROW( runtime_error, "Error: quality data for read " << m_readId << " (2) does not correspond to sequence length" );
            m_readData[0] += qf;
            m_readData[1] += qr;
        }
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////
// CSolexaFileReader

CSolexaFileReader::CSolexaFileReader( const string& readIdFile, const string& readDataFile1, const string& readDataFile2 )
{
    SetStream( m_idStream, readIdFile );
    SetStream( m_dataStream1, readDataFile1 );
    if( readDataFile2.length() )
        SetStream( m_dataStream2, readDataFile2 );
}

bool CSolexaFileReader::NextRead() 
{
    m_readId = m_readData[0] = m_readData[1] = "";
    string buff;
    while( getline( m_idStream->SetStream(), buff ) ) {
        buff = Trim( buff );
        if( buff.length() == 0 || buff[0] == '#' ) continue;
        istringstream in( buff );
        in >> m_readId;
        if( ! getline( m_dataStream1->SetStream(), m_readData[0] ) ) 
            THROW( runtime_error, "Failed to read solexa data from " << m_dataStream1->GetName() << " for read " << m_readId << ": inconsistent input?" );
        if( !x_PairedReads() ) return true;
            if( ! getline( m_dataStream2->SetStream(), m_readData[1] ) ) 
            THROW( runtime_error, "Failed to read solexa data from " << m_dataStream2->GetName() << " for read " << m_readId << ": inconsistent input?" );
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////
// CFastqFileReader

CFastqFileReader::CFastqFileReader( bool quality, bool colorspace, const string& readDataFile1, const string& readDataFile2 )
{
    if( colorspace ) {
        if( quality ) THROW( logic_error, "Oops... don't know how to use quality scores in colorspace" );
        m_coding = CSeqCoding::eCoding_colorsp;
    } else m_coding = quality ? CSeqCoding::eCoding_ncbiqna : CSeqCoding::eCoding_ncbi8na;

    SetStream( m_dataStream1, readDataFile1 );
    if( readDataFile2.length() )
        SetStream( m_dataStream2, readDataFile2 );
}

bool CFastqFileReader::x_FetchIdLine( TFileStream& stream, string& dest, char type, const string& compare, bool allowComponentDifference )
{
    string buff;
    while( getline( stream->SetStream(), buff ) ) {
        if( IgnoreHashLines() && buff.length() && buff[0] == '#' ) continue;
        x_TrimTrailingSpaces( buff );
        if( buff.length() == 0 ) continue; // it's better to skip empty lines to avoid complications with end of file
        if( buff[0] != type ) 
            THROW( runtime_error, "Error: bad input line of [" << buff << "] in " << stream->GetName() << ": expected to start with ``" << type << "''" );
        const char * e = buff.c_str();
        while( *e && !isspace( *e ) ) ++e;
        dest.assign( buff.c_str() + 1, e );
        if( compare.length() && dest.length() && compare != dest ) {
            if( allowComponentDifference && compare.length() > 2 && dest.length() == compare.length() ) {
                int l = compare.length() - 2;
                if( compare.substr( 0, l+1 ) == dest.substr( 0, l+1 ) && 
                        compare[l+1] == '1' && dest[l+1] == '2' && strchr( "/_.", compare[l] ) )
                            return m_clipReadId = true;
            }
            THROW( runtime_error, "Error: bad read id in file " << stream->GetName() << ": expected " << compare << ", got " << dest );
        }
        return true;
    }
    return false;
}

bool CFastqFileReader::x_FetchReadData( TFileStream& stream, const string& id, string& dest, int expectedLen )
{
    while( true ) {
        int c = stream->SetStream().get();
        if( stream->SetStream().fail() ) return false;
        stream->SetStream().putback( c );
        switch( m_coding ) {
        case CSeqCoding::eCoding_colorsp:
        case CSeqCoding::eCoding_ncbi8na:
            if( c == '>' ) return true;
            else {
                string buff;
                do getline( stream->SetStream(), buff );
                while( IgnoreHashLines() && buff.length() && buff[0] == '#' );
                x_TrimTrailingSpaces( buff );
                dest += buff;
            }
            break;
        case CSeqCoding::eCoding_ncbiqna:
            if( expectedLen == 0 ) {
                if( c == '+' ) {
                    if( dest.length() == 0 ) THROW( runtime_error, "Error: No sequence data in fastq file " << stream->GetName() << " for read " << id );
                    string foo;
                    if( !x_FetchIdLine( stream, foo, '+', id, false ) ) 
                        THROW( runtime_error, "Error: Unexpected end of file " << stream->GetName() << ": no quality data for read " << id );
                    string qual;
                    x_FetchReadData( stream, id, qual, dest.length() );
                    ASSERT( qual.length() == dest.length() );
                    dest += qual;
                    return true;
                } else {
                    string buff;
                    getline( stream->SetStream(), buff );
                    x_TrimTrailingSpaces( buff );
                    dest += buff;
                }
            } else {
                do {
                    string buff;
                    getline( stream->SetStream(), buff );
                    x_TrimTrailingSpaces( buff );
                    dest += buff;
                } while( expectedLen > (int)dest.length() );
                if( expectedLen < (int)dest.length() ) 
                    THROW( runtime_error, "Error: Quality data [" << dest << "] are longer then expected (" << expectedLen << ") in file " << stream->GetName() << " read " << id );
                if( expectedLen == (int)dest.length() ) return true;
            }
            break;
        default: THROW( logic_error, "Internal Error: Wrong coding in " __FILE__ ": " << __LINE__ );
        }
    }
    return false;
}

bool CFastqFileReader::NextRead()
{
    m_readId = m_readData[0] = m_readData[1] = "";
    m_clipReadId = false;
    char type = m_coding == CSeqCoding::eCoding_ncbiqna ? '@' : '>';
    while( x_FetchIdLine( m_dataStream1, m_readId, type, "" ) ) {
        x_FetchReadData( m_dataStream1, m_readId, m_readData[0], 0 );
        if( x_PairedReads() ) {
            string foo;
            if( !x_FetchIdLine( m_dataStream2, foo, type, m_readId, true ) )
                THROW( runtime_error, "Error: failed to fetch data for " << m_readId << " from " << m_dataStream2->GetName() << ": truncated file?" );
            x_FetchReadData( m_dataStream2, foo, m_readData[1], 0 );
            if( m_clipReadId ) m_readId.resize( m_readId.length() - 2 );
        }
        if( m_coding == CSeqCoding::eCoding_colorsp ) {
            for( int i = 0; i < 2; ++i ) {
                size_t dot = m_readData[i].find('.');
                if( dot != string::npos ) {
                    cerr << "* Warning: read " << m_readId << " component " << i << " is clipped to " << dot << " bases\n";
                    m_readData[i].resize( dot );
                }
            }
        }
        return true;
    }
    return false;
}

