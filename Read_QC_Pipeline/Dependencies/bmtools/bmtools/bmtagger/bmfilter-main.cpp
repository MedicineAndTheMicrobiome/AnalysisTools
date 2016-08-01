#include "bmtagger-version.hpp"
#include "cprogressindicator.hpp"
#include "csrafastqreader.hpp"
//#include "cbitmaskaccess.hpp"
#include "cshortreader.hpp"
#include "creadtagger.hpp"
#include "csamrecords.hpp"
#include "cgetopt.hpp"

#include <limits>

USING_OLIGOFAR_SCOPES;

int main( int argc, char ** argv )
{
#if WITH_SRA
    set<string> accession;
    bool spotIdOnly = false;
    bool sra2fasta = false;
    bool sraFilterBad = false;
    bool sraFilterCond = false;
    bool sraFilterHuman = false;
    bool sraTagFiltered = false;
    Uint4 sraFlags = 0;
    set<string> sraSpotIdS;
    set<int> sraSpotIds;
#endif
    set<string> bmask;
    list<string> in1;
    list<string> in2;
    bool mmap = false;

    CReadTagger tagger;
    tagger.SetShortSeqSwitch( numeric_limits<int>::max()/2 );
    tagger.SetChopLength( 32 );

    CGetOpt getopt( argc, argv );
    getopt
        .AddGroup( "general" )
        .AddHelp( argv[0], "Filter reads based on word bitmask matches" )
        .AddVersion( argv[0], BMTAGGER_VERSION, __DATE__ )
        .AddGroup( "input" )
        .AddArg( tagger.SetQualityChannels(), 'q', "quality-channels", "Number of quality channers for reads (0|1)" )
        .AddList( in1, '1', "read-1", "Fasta or fastq (for -q1) file with reads, may be repeated" )
        .AddList( in2, '2', "read-2", "Fasta or fastq (for -q1) file with read pair mates, if used should be repeated as many times as -1 is" )
        .AddSet( bmask, 'b', "word-bitmask", "Word bitmask file (may be repeated)" )
        .AddFlag( mmap, 'M', "use-mmap", "Use mmap for word bitmask (slow unless used for few reads; intended for debug)" )
#if WITH_SRA
        .AddSet( accession, 'A', "accession", "SRA accession to process (may be repeated)" )
        .AddSet( sraSpotIdS, 'i', "spot-id", "SRA spot-ids to process (may be repeated)" )
#endif
        .AddGroup( "input-filters" )
        .AddArg( tagger.SetMaxAmbiguities(), 'a', "max-ambiguities", "Maximal number of ambiguities per word" )
        .AddFlag( tagger.SetClipLowercase(), 'l', "clip-lowercase", "Should lowercase head and tail of each read be clipped" )
        .AddArg( tagger.SetClipNwindow(), 'N', "clip-N-win", "Clip sequence head or tail as long as it has at least one N per this long window" )
        .AddArg( tagger.SetClipQuality(), 'Q', "clip-quality", "Clip sequence head or tail with quality lower then this (for fastq input)" )
#if WITH_SRA
        .AddFlag( sraFilterBad, 0, "sra-filter-bad", "Don't process bad reads" )
        .AddFlag( sraFilterCond, 0, "sra-filter-cond", "Don't process conditional reads" )
        .AddFlag( sraFilterHuman, 0, "sra-filter-human", "Don't process human reads" )
#endif
        .AddGroup( "output" )
        .AddArg( tagger.SetOutputBasename(), 'o', "output", "Output base name (suffixes will be added to it)" )
#if WITH_SRA
        .AddFlag( spotIdOnly, 'I', "spot-id-only", "Don't use `accession.' for output (only available with single SRA accession)" )
        .AddFlag( sraTagFiltered, 0, "sra-tag-filtered", "Tag SRA-filtered hits as human" )
        .AddFlag( sra2fasta, 0, "dump-as-fasta", "Write input reads to output fasta .fa file" )
        .AddFlag( (Uint4)CSRAFastqReader::fReadDb_useFastqReader, sraFlags,  0, "sra-use-fastq-reader", "Use FastqReader from SRA toolkit" )
        .AddFlag( (Uint4)CSRAFastqReader::fReadDb_readClipInfo,   sraFlags,  0, "sra-read-clipinfo", "Read clip info columns (not with FastqReader)" )
#endif
        .AddFlag( (Uint8)CReadTagger::fAction_tag, tagger.SetActions(), 'T', "tag", "Produce .tag file" )
        .AddFlag( (Uint8)CReadTagger::fAction_post, tagger.SetActions(), 'P', "post", "Produce .short?.fa and .long?.fa files" )
        .AddFlag( (Uint8)CReadTagger::fAction_reportNoseq, tagger.SetActions(), 'R', "report", "Produce .report file" )
        .AddFlag( (Uint8)CReadTagger::fAction_postClipped, tagger.SetActions(), 'z', "post-clipped", "Put clipped versions of sequences to output .fa files" )
        .AddArg( tagger.SetOutComplexityFilterCutoff(), 'F', "complexity", "Set complexity filter cutoff" )
        .AddArg( tagger.SetShortSeqSwitch(), 'L', "short-seq", "Set sequence length to consider it as short for postprocessing" )
        .AddArg( tagger.SetNoPostSeqLength(), 'n', "no-post-len", "Set longest sequence length to ignore postprocessing" )
        .AddArg( tagger.SetChopLength(), 'c', "chop-length", "Set length to chop short sequences to" )
        .AddArg( tagger.SetChopStep(), 's', "chop-step", "Set step by which to chop short sequences" )
        .AddArg( tagger.SetMaskComplexityEarly(), 'm', "mask-early", "Set mask low complexity before applying heuristics" )
        .AddArg( tagger.SetPostLowComplexity(), 'Z', "post-low-complexity", "Should 'unknown' low complexity reads be sent to post processing" )
        .AddGroup( "heuristics" )
        .AddArg( tagger.SetMinWordsCutoff(), 0, "heur-min-words", "Set minimal word count to apply heuristics" )
        .AddArg( tagger.SetLongWordsSwitch(), 0, "heur-many-words", "Set number of good words which switches watermarks (long/short)" )
        .AddArg( tagger.SetHeurMatchCountLongPct(), 0, "heur-count-long-pct", "Set watermarks for matched word count for long sequences, int % of good words" )
        .AddArg( tagger.SetHeurMatchCountShortPct(), 0, "heur-count-short-pct", "Set watermarks for matched word count for short sequences, int % of good words" )
        .AddArg( tagger.SetHeurMatchLongestLongPct(), 0, "heur-run-long-pct", "Set watermarks for longest match run for long sequences, int % of good words" )
        .AddArg( tagger.SetHeurMatchLongestShortPct(), 0, "heur-run-short-pct", "Set watermarks for longest match run for short sequences, int % of good words" )
        .AddArg( tagger.SetNegligibleSeqLength(), 0, "heur-negligible-length", "Set cutoff for sequences to consider - these and shorter (after clipping) will be marked as foreign" ) 
        .Parse();
    if( getopt.Done() ) return getopt.GetResultCode();
    
    if( in2.size() ) {
        if( in2.size() != in1.size() ) 
            THROW( runtime_error, "if -2 is used it should be repeated exactly as many times as -1 is, usage is (" << in1.size() << ", " << in2.size() << ")" );
    }
    
    if( tagger.GetOutputBasename().length() == 0 )
        THROW( runtime_error, "USER error: no output base filename is given: option -o is required" );

    if( bmask.size() == 0 ) 
        THROW( runtime_error, "USER error: no bitmask files are provided; please use -b option" );

    ITERATE( set<string>, bma, bmask ) {
        cerr << "* Attaching " << *bma << "...";
        cerr.flush();
        unsigned flags = mmap ?CBitmaskAccess::fUseMmap : 0;
        auto_ptr<CBitmaskAccess> b( new CBitmaskAccess() );
        b->AttachFile( *bma, flags ); // mmap ? CBitmaskAccess::eOpen_mmap : CBitmaskAccess::eOpen_read ) );
        cerr << "\b\b\b (pattern = "
            << "0b" << CBitHacks::AsBits( b->GetBasePattern(), b->GetWindowLength() ) 
            << " of len " << b->GetWindowLength() << " using " << (2*b->GetWordBases()) << " bits)\n";
        tagger.AddBitmaskAccess( b.release(), true );
    }

#if WITH_SRA
    if( spotIdOnly && accession.size() > 1 )
        THROW( runtime_error, "USER Error: it is confusing to use -I with multiple -A" );
    auto_ptr<ofstream> ofa(0);
    if( sra2fasta ) ofa.reset( new ofstream( (tagger.GetOutputBasename() + ".fa").c_str() ) );
    ITERATE( set<string>, acc, accession ) {
        CSRAFastqReader sraReader( *acc, sraFlags );
        sraReader.AssignFlags( sraFlags );
        sraReader.SetFlags( CSRAFastqReader::fReadId_includeAccession, !spotIdOnly );
        sraReader.SetRqQuality( tagger.GetQualityChannels() );
        sraReader.Rewind();
        ITERATE(set<string>, x, sraSpotIdS) sraSpotIds.insert( FromStr<int>( *x ) );
        sraReader.SetSpotIdList( sraSpotIds.begin(), sraSpotIds.end() );
        CProgressIndicator p( "Processing SRA " + *acc );
        for( ; sraReader.Exists() && sraReader.FetchRead() ; sraReader.NextRead() ) {
            if( sra2fasta ) {
                for( int i = 0; i < sraReader.GetComponents(); ++i ) {
                    static const char * flt[] = { "pass", "reject", "criteria", "redacted" };
                    *ofa << ">" << sraReader.GetReadId() << "/" << (i+1) 
                        << " len=" << sraReader.GetRead( i ).length()
                        << " clip={" << (sraReader.GetClip(i).first + 1) << ".." << sraReader.GetClip(i).second << "}"
                        << " " << flt[sraReader.GetReadFilterValue()] << "\n"; 
                    *ofa << sraReader.GetRead( i ) << "\n";
                }
            }
            bool filtered = false;
            switch( sraReader.GetReadFilterValue() ) {
                case CSRAFastqReader::eReadFilter_pass: break;
                case CSRAFastqReader::eReadFilter_reject: if( sraFilterBad ) filtered = true; break;
                case CSRAFastqReader::eReadFilter_criteria: if( sraFilterCond ) filtered = true; break;
                case CSRAFastqReader::eReadFilter_redacted: if( sraFilterHuman ) filtered = true; break;
                default: THROW( runtime_error, "Unexpected value of sraReader.GetReadFilterValue(): " << sraReader.GetReadFilterValue() );
            }
            if( filtered ) { if( sraTagFiltered ) tagger.ForceHuman( filtered ); else continue; }
            if( sraReader.GetComponents() == 0 ) continue;
            tagger.SetComponents( sraReader.GetComponents() );
            for( int i = 0; i < sraReader.GetComponents(); ++i ) {
                const CSRAFastqReader::TRange& range = sraReader.GetClip( i );
                //cout << ">" << sraReader.GetReadId() << "." << i << "\t" << sraReader.GetRead(i).length() << "\t" << range.GetFrom() << ".." << range.GetTo() << "\t" << tagger.GetTagValue() << "\n";
                //cout << sraReader.GetRead( i ) << "\n";
                if( tagger.GetQualityChannels() ) {
                    ASSERT( sraReader.GetQual( i ).length() == sraReader.GetRead( i ).length() );
                    string rdata = sraReader.GetRead( i ) + sraReader.GetQual( i );
                    tagger.ProcessReadData( sraReader.GetReadId().c_str(), rdata.c_str(), 0/*range.first*/, range.second == range.first ? 0 : range.second );
                } else {
                    tagger.ProcessReadData( sraReader.GetReadId().c_str(), sraReader.GetRead( i ).c_str(), 0/*range.first*/, range.second == range.first ? 0 : range.second );
                }
            }
            p.Increment();
        }
        p.Summary();
    }
#endif

    CProgressIndicator p( "Processing files" );
    list<string>::const_iterator i1 = in1.begin();
    list<string>::const_iterator i2 = in2.begin();
    while( i1 != in1.end() ) {
        if( i1 != in1.begin() ) cerr << "\n"; //p.Summary();
        CReaderFactory rfactory;
        rfactory.SetQualityChannels( tagger.GetQualityChannels() );
        rfactory.SetReadDataFile1( *i1 );
        if( i2 != in2.end() ) { rfactory.SetReadDataFile2( *i2 ); tagger.SetComponents( 2 ); }
        else tagger.SetComponents( 1 );

        auto_ptr<IShortReader> reader( rfactory.CreateReader() );
        if( reader.get() == 0 ) {
            THROW( runtime_error, "Failed to create short read file reader" );
            return 10;
        }

        while( reader->NextRead() ) {
            tagger.ProcessReadData( reader->GetReadId(), reader->GetReadData(0) );
            if( tagger.GetComponents() > 1 )
                tagger.ProcessReadData( reader->GetReadId(), reader->GetReadData(1) );
            p.Increment();
        }
        ++i1;
        if( i2 != in2.end() ) ++i2;
    }
    tagger.PurgeRead();
    tagger.Complete();
    p.Summary();

    return 0;
}

