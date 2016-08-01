#include "cgetopt.hpp"
#include "cseqtobitmask.hpp"
#include "cseqwcbitmask.hpp"
#include "bmtagger-version.hpp"
#include "bmask-tmpl.hpp"

USING_OLIGOFAR_SCOPES;


int main( int argc, char ** argv )
{
    string seqdb;
    string wmask;
    string gilist;
    string featfile;
    string input;
    int wsize = 18;
    int wstep = 1;
    int maxamb = 0;
    unsigned pattern = 0;
    bool fastaParseIDs = false;
    bool quiet = false;

    bool test = false;
    bool mmap = false;
    bool diff = false;
    bool slow = false;
    bool compress = false;
    bool extraCompress = false;
    int minor = 2;
    int packPrefixBits = 26;
    int packOffsetBits = 34;
    int packCountBits = 10;
    int maxWordCount = 0;

    CGetOpt getopt( argc, argv );
    getopt
        .AddGroup( "general" )
        .AddHelp( argv[0], "Build word mask for subject database" )
        .AddVersion( argv[0], BMTAGGER_VERSION, __DATE__ )
        .AddFlag( quiet, 'q', "quiet", "Do not show progress indicators" )
        .AddGroup( "files" )
        .AddArg( seqdb, 'd', "fasta-file", "Input fasta or blastdb file" )
        .AddArg( wmask, 'o', "output-file", "Output word bitmask file" )
        .AddArg( gilist, 'l', "gi-list", "Set gi list for blastdb file" )
        .AddArg( input, 'i', "input-file", "Set word bitmask file as input" )
        .AddArg( fastaParseIDs, 0, "fasta-parse-ids", "Parse FASTA ids (becomes broken if ranges are used)" )
        .AddGroup( "hashing parameters" )
        .AddArg( wsize, 'w', "word-size", "Word size to use" )
        .AddArg( wstep, 'S', "word-step", "Step (stride size) to use" )
        .AddArg( maxamb, 'A', "max-amb", "Maximal number of ambiguities to count" )
        .AddArg( pattern, 'p', "pattern", "Set pattern to use with discontiguous words, 0x or 0b prefix may be used for hex or bin (-w## will be ignored)" )
        .AddArg( maxWordCount, 'W', "max-word-count", "Do not include words counted more than this value (for 16-mers or less)" )
        .AddGroup( "output" )
        .AddArg( minor, 'v', "version", "Create this version of file (0-2)" )
        .AddFlag( compress, 'z', "compress", "Compress bitmask (requires version 2)" )
        .AddFlag( extraCompress, 'Z', "extra-compress", "Compress bitmask (requires version 2) looking for duplicate extension sets" )
        .AddArg( packPrefixBits, 0, "pack-prefix-bits", "Bits to use for compression prefix" )
        .AddArg( packOffsetBits, 0, "pack-offset-bits", "Number of bits in table to use for data segment offset" )
        .AddArg( packCountBits, 0, "pack-count-bits", "Number of bits to reserve for entry count within segment" )
        .AddGroup( "other" )
        .AddFlag( mmap, 0, "mmap", "Memory map source file instead of reading" )
        .AddFlag( diff, 0, "diff", "Diff source and result before writing, repport differences" )
        .AddFlag( slow, 0, "slow", "Slow copy (using query API - to check query api" )
        .AddFlag( test, 0, "bit-test", "Test bit operations" )
        .Parse();
    if( getopt.Done() ) return getopt.GetResultCode();

#define CHECK( ii, jj, x ) \
    if( CBitHacks::GetBits( &data[0], jj, ii ) != x ) \
        THROW( logic_error, "Expected from GetBits( ptr, " << jj << ", " << ii << ") value 0x" << hex << x << ", got 0x" << CBitHacks::GetBits( &data[0], jj, ii ) << " for i = " << i << ", j = " << j << ", and w = 0x" << hex << w )
#define PRINT \
        cout << setw(3) << i << ", " << setw(3) << j << setw(0) << ": "; ITERATE( vector<unsigned char>, d, data ) cout << CBitHacks::AsBits( *d, 8 ); cout << "\n"

    if( extraCompress ) compress = true;

    if( test ) {
        const int bits = 128;
        for( int i = 1; i < 63; ++i ) {
            Uint8 w = CBitHacks::WordFootprint<Uint8>(i);
            for( int j = 0; j <= bits - i; ++j ) {
                vector<unsigned char> data(bits/8);
                fill( data.begin(), data.end(), 0 );
                CBitHacks::SetBits( &data[0], j, i, ~Uint8(0) );
                //PRINT;
                ASSERT( CBitHacks::GetBits( &data[0], j, i ) == w );
                if( j ) ASSERT( CBitHacks::GetBits( &data[0], j-1, i+1 ) == w );
                if( j + i < bits ) ASSERT( CBitHacks::GetBits( &data[0], j, i+1 ) == (w << 1) );
            }
        }
        for( int i = 1; i < 64; ++i ) {
            for( int j = 0; j <= bits - i; ++j ) {
                vector<unsigned char> data(bits/8);
                fill( data.begin(), data.end(), -1 );
                CBitHacks::SetBits( &data[0], j, i, 0 );
                PRINT;
            }
        }
        return 0;
    }

    if( seqdb.length() == 0 || seqdb == "-" ) {
        if( wmask.length() == 0 && input.length() == 0 ) 
            THROW( runtime_error, "Either input or output file name is required" );
        seqdb = "/dev/stdin";
    } else if( wmask.length() == 0 ) {
        wmask = seqdb + ".wbm";
    }

    if( wstep != 1 ) THROW( runtime_error, "Sorry, only wstep=1 is supported now" );
    
    if( minor < 2 && compress) {
        compress = false;
        cerr << "WARNING: compression is ignored in compatibility mode\n";
    }

    if( input.length() ) {
        if( seqdb.size() )
            cerr << "WARNING: can't use seqdb when repacking data - ignoring file(s) " << seqdb << "*\n";
        cerr << "INFO: attaching file " << input << "\n";
        Uint4 flags = mmap ? CBitmaskAccess::fUseMmap|CBitmaskAccess::fMmapSequential : 0;

        CBitmaskAccess bmaccess;
        CProgressIndicator pr( "Reading file "+input+", Mb" );
        bmaccess.AttachFile( input, flags, quiet?0:&pr );
        if( !quiet ) pr.Summary();

        CBitmaskWriter writer;
        writer.SetFileMinor(minor);
        writer.SetMaxAmb( bmaccess.GetMaxAmb() );
        writer.SetPattern( bmaccess.GetBasePattern() );
        if( compress ) {
            CProgressIndicator p( "Repacking compressed bitmask" );
            if( minor < 2 ) { 
                minor = 2;
                cerr << "Warning: using file version 0." << minor << ".0 for compression";
            } 
            CVarBmPackConfig cfg( packPrefixBits, bmaccess.GetWordBits() - packPrefixBits, packOffsetBits, packCountBits );
            typedef CPackedBitmaskDataAppender<CVarBmPackConfig> TAppender;
            TAppender bmappender( cfg );
            bmappender.SetUseHashOffsets( extraCompress );
            if( slow ) {
                for( Uint8 w = 0; w <= bmaccess.GetWordMask(); ++w ) {
                    if( bmaccess.HasExactWord( w ) )
                        bmappender.AddExactWord( w );
                    if( !quiet ) p.Increment();
                }
            } else {
                CBitmaskAccess::C_Builder<TAppender> builder( bmappender, quiet?0:&p );
                bmaccess.ForEachSet( builder );
            }
            bmappender.Complete();
            if( !quiet ) p.Summary();
            if( diff ) {
                CProgressIndicator pd( "Running diff" );
                CPackedBitmaskDataReader<CVarBmPackConfig> bmreader( bmappender.GetConfig(), bmappender.GetTable(), bmappender.GetData() );
                for( Uint8 w = 0; w <= bmaccess.GetWordMask(); ++w ) {
                    bool a = bmaccess.HasExactWord( w );
                    bool b = bmreader.HasExactWord( w );
                    if( a != b ) {
                        //DEBUG( DISPLAY( a ) << DISPLAY( b ) << DISPLAYh( w ) << DISPLAY( w ) << bmreader.GetSuffixListCount( w ) );
                        cerr << bmreader.GetConfig().Word2Suffix( w ) << ":\t";
                        bmreader.GetSuffixList( w, ostream_iterator<Uint8>( cerr, ", " ) );
                        cerr << "\n";
                    }
                    pd.Increment();
                }
                if( !quiet ) pd.Summary();
            }
            //DEBUG( DISPLAY( bmappender.GetConfig().GetTableBytesPacked() ) << DISPLAY( bmappender.GetDataSize() ) << DISPLAY( bmappender.GetBitCount() ) << DISPLAY( bmappender.GetAddCount() ) << DISPLAY( bmappender.GetSavedBytes() ) );
            writer.Write( wmask, bmappender, bmappender.GetDataSize() );
        } else {
            CProgressIndicator p( "Repacking flat bitmask" );
            typedef CFlatBitmaskData<Uint4> TAppender;
            TAppender bmappender( CFlatBitmaskData<Uint4>::Allocate( bmaccess.GetWordBases() ), bmaccess.GetWordBases() );
            if( slow ) {
                for( Uint8 w = 0; w <= bmaccess.GetWordMask(); ++w ) {
                    if( bmaccess.HasExactWord( w ) )
                        bmappender.AddExactWord( w );
                    if( !quiet ) p.Increment();
                }
            } else {
                CBitmaskAccess::C_Builder<TAppender> builder( bmappender, quiet?0:&p );
                bmaccess.ForEachSet( builder );
            }
            if( !quiet ) p.Summary();
            CProgressIndicator pw( "Writing flat bitmask to " + wmask );
            writer.Write( wmask, bmappender, quiet?0:&pw );
            if( !quiet ) pw.Summary();
        }
    } else {
        if( pattern == 0 ) pattern = CBitHacks::WordFootprint<Uint8>( wsize );
        //DEBUG( DISPLAY( pattern ) << DISPLAY( wsize ) );

        typedef CFlatBitmaskData<Uint4> TAppender;
        CBitmaskWriter writer;
        writer.SetPrefixBits( 0 );
        writer.SetFileMinor(minor);
        writer.SetMaxAmb( maxamb );
        writer.SetPattern( pattern );

        TAppender bmappender( CFlatBitmaskData<Uint4>::Allocate( writer.GetWordBases() ), writer.GetWordBases() );
        auto_ptr<CSeqVecProcessor::ICallback> builder(0);
        if( maxWordCount > 0 ) 
            builder.reset( new CSeqWcBitmask<TAppender>( &bmappender, pattern, maxamb, quiet ) );
        else
            builder.reset( new CSeqToBitmask<TAppender>( &bmappender, pattern, maxamb, quiet ) );
        //CSeqToBitmask<TAppender> builder( &bmappender, pattern, maxamb );
        cerr << "* Info: using 4^" << writer.GetWordBases() << " bits for ";
        if( writer.IsDiscontiguous() )
            cerr << "discontiguous words " << CBitHacks::AsBits( writer.GetBasePattern() );
        cerr << " of " << writer.GetWindowLength() << " bases\n";
        if( maxWordCount > 0 ) 
            cerr << "* Info: ignoring overrepresented words (for counts above " << maxWordCount << ")\n";

        CSeqVecProcessor processor;
        processor.AddCallback( 0, builder.get() );
#ifndef STANDALONE
        if( gilist.size() ) processor.SetGiListFile( gilist );
        processor.SetFastaReaderParseId( fastaParseIDs );
#endif
        processor.Process( seqdb );
        if( CSeqWcBitmask<TAppender> * x = dynamic_cast<CSeqWcBitmask<TAppender>*>( builder.get() ) ) 
            x->CommitToBitmask( maxWordCount );

        builder.reset(0);

        CProgressIndicator pw( "Writing " + string( compress ? "compressed" : "flat" ) + " bitmask to " + wmask );
        if( compress ) {
            typedef CPackedBitmaskDataAppender<CVarBmPackConfig> TPacker;
            CVarBmPackConfig cfg( packPrefixBits, writer.GetWordBits() - packPrefixBits, packOffsetBits, packCountBits );
            TPacker bmpacker( cfg );
            bmpacker.SetUseHashOffsets( extraCompress );
            CProgressIndicator p( "Compressing bitmask" );
            CBitmaskAccess::C_Builder<TPacker> packer( bmpacker, quiet?0:&p );
            bmappender.ForEachSet( packer );
            bmpacker.Complete();
            if( !quiet ) p.Summary();
            //DEBUG( DISPLAY( bmpacker.GetConfig().GetTableBytesPacked() ) << DISPLAY( bmpacker.GetDataSize() ) << DISPLAY( bmpacker.GetBitCount() ) << DISPLAY( bmpacker.GetAddCount() ) << DISPLAY( bmpacker.GetSavedBytes() ) );
            writer.Write( wmask, bmpacker, bmpacker.GetDataSize(), quiet?0:&pw );
        } else {
            writer.Write( wmask, bmappender, quiet?0:&pw );
        }
        if( !quiet ) pw.Summary();
    }

    return 0;
}

