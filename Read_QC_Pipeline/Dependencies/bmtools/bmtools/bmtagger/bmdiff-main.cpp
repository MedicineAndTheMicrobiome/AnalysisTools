#include "cgetopt.hpp"
#include "cprogressindicator.hpp"
#include "bmask-tmpl.hpp"

USING_OLIGOFAR_SCOPES;

void PrintBitmaskInfo( ostream& o, const CBitmaskAccess& bm, const string& fname )
{
    o 
        << "filename:  " << fname << "\n"
        << "wordsize:  " << bm.GetWordBases() << "\n"
        << "wordstep:  " << bm.GetWinStep() << "\n"
        << "winlength: " << bm.GetWindowLength() << "\n"
        << "pattern:   0b" << CBitHacks::AsBits( bm.GetBasePattern(), bm.GetWindowLength() ) << "\n"
        ;
}

int main( int argc, char ** argv )
{
    bool mmap = false;
    bool force = false;

    CGetOpt getopt( argc, argv );
    getopt
        .AddGroup( "general" )
        .AddHelp( "usage: bmdiff <bitmask1> <bitmask2>", "Compare two bitmasks created by bmtool, print statistics" )
        .AddVersion( argv[0], "0.0.0", __DATE__ )
        .AddGroup( "options" )
        .AddArg( mmap, 'm', "mmap", "Memory map input files" )
        .AddFlag( force, 'f', "force", "Force comparison if patterns or window steps differ" )
        .Parse();
    if( getopt.Done() ) return getopt.GetResultCode();

    if( getopt.GetArgCount() != getopt.GetArgIndex() + 2 )
        THROW( runtime_error, "Need exactly two positional arguments" );

    Uint4 flags = mmap ? CBitmaskAccess::fUseMmap|CBitmaskAccess::fMmapSequential : 0;

    string file1(getopt.GetArg( getopt.GetArgIndex() + 0 ));
    string file2(getopt.GetArg( getopt.GetArgIndex() + 1 ));
    CBitmaskAccess bm1, bm2;
    CProgressIndicator pr1( "Attaching " + file1, ", Mb" );
    bm1.AttachFile( file1, flags, &pr1 );
    pr1.Summary();
    CProgressIndicator pr2( "Attaching " + file2, ", Mb" );
    bm2.AttachFile( file2, flags, &pr2 );
    pr2.Summary();

    if( bm1.GetWordBases() != bm2.GetWordBases() ) {
        THROW( runtime_error, "Can't compare bitmasks with different sizes: " << bm1.GetWordBases() << " and " << bm2.GetWordBases() );
    }
    if( bm1.GetBasePattern() != bm2.GetBasePattern() ) {
        if( force ) 
            cerr << "Warning: patterns of bitmasks differ: " << hex << bm1.GetBasePattern() << " and " << bm2.GetBasePattern() << dec << "\n";
        else 
            THROW( runtime_error, "Can't compare bitmasks with different pattern: "  << hex << bm1.GetBasePattern() << " and " << bm2.GetBasePattern() << dec );
    }
    if( bm1.GetWinStep() != bm2.GetWinStep() ) {
        if( force ) 
            cerr << "Warning: word steps of bitmasks differ: " << bm1.GetWinStep() << " and " << bm2.GetWinStep() << "\n";
        else 
            THROW( runtime_error, "Can't compare bitmasks with different word steps: "  << hex << bm1.GetWinStep() << " and " << bm2.GetWinStep() );
    }

    Uint8 stat[4];
    fill( stat, stat+4, 0 );
    CProgressIndicator pd( "Running diff", "megabits" );
    const Uint8 step = 1024*1024;
    const Uint8 mask = step-1;
    pd.SetFinalValue( bm1.GetWordMask()/step );
    for( Uint8 w = 0; w <= bm1.GetWordMask(); ++w ) {
        bool a = bm1.HasExactWord( w );
        bool b = bm2.HasExactWord( w );
        int x = int(a) | (int(b)<<1);
        stat[x] ++;
        if( (w&mask) == mask ) pd.Increment();
    }
    pd.Summary();

    cout << "Compared bitmasks\n"
         << "=================" ;
    cout << "\nFirst bitmask:\n";
    PrintBitmaskInfo( cout, bm1, file1 );
    cout << "\nSecond bitmask:\n";
    PrintBitmaskInfo( cout, bm2, file2 );
    Uint8 sum = stat[0] + stat[1] + stat[2] + stat[3];
    cout << "\n--------------------------------------------------\n";
    cout << "Words absent in both inputs:  " << setw(12) << right << stat[0] << " " << setw(6) << fixed << setprecision(2) << (100.0*stat[0])/sum << "%\n";
    cout << "Words present in first only:  " << setw(12) << right << stat[1] << " " << setw(6) << fixed << setprecision(2) << (100.0*stat[1])/sum << "%\n";
    cout << "Words present in second only: " << setw(12) << right << stat[2] << " " << setw(6) << fixed << setprecision(2) << (100.0*stat[2])/sum << "%\n";
    cout << "Words present in both inputs: " << setw(12) << right << stat[3] << " " << setw(6) << fixed << setprecision(2) << (100.0*stat[3])/sum << "%\n";
    cout << "Total words:                  " << setw(12) << right << sum << "\n";
    return 0;
}

