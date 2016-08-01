#include "cgetopt.hpp"
#include "cseqtobitmask.hpp"
#include "cseqwcbitmask.hpp"
#include "bmtagger-version.hpp"
#include "bmask-tmpl.hpp"

USING_OLIGOFAR_SCOPES;


int main( int argc, char ** argv )
{
    string input;
    string output;
    bool mmap = false;

    CGetOpt getopt( argc, argv );
    getopt
        .AddGroup( "general" )
        .AddHelp( argv[0], "Dump word mask as UIPACna" )
        .AddVersion( argv[0], BMTAGGER_VERSION, __DATE__ )
        .AddGroup( "files" )
        .AddArg( output, 'o', "output-file", "Output words file" )
        .AddArg( input, 'i', "input-file", "Set word bitmask file as input" )
        .AddGroup( "other" )
        .AddFlag( mmap, 0, "mmap", "Memory map source file instead of reading" )
        .Parse();
    if( getopt.Done() ) return getopt.GetResultCode();
    if( input.size() == 0 ) {
        cerr << "Input file is required\n";
        return 64;
    }
    ostream * out = &cout;
    ofstream fout;
    if( output.size() ) {
        fout.open( output.c_str() );
        out = &fout;
    }

    cerr << "INFO: attaching file " << input << "\n";
    Uint4 flags = mmap ? CBitmaskAccess::fUseMmap|CBitmaskAccess::fMmapSequential : 0;

    CBitmaskAccess bmaccess;
    CProgressIndicator pr( "Reading file "+input+", Mb" );
    bmaccess.AttachFile( input, flags, &pr );
    pr.Summary();

    for( Uint8 w = 0; w <= bmaccess.GetWordMask(); ++w ) {
        if( bmaccess.HasExactWord( w ) ) {
            string x( bmaccess.GetWordBases(), '.' );
            for( Uint8 p = bmaccess.GetWordBases(), ww = w, i = 0; p > 0; ++i, --p, (ww>>=2) ) 
                x[i]= "ACGT"[ww&3];
            *out << x << "\n";
        }
    }
    return 0;
}

