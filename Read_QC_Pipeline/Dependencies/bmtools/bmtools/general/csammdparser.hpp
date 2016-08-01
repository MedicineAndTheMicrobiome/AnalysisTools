#ifndef OLIGOFAR_CSAMMDPARSER__HPP
#define OLIGOFAR_CSAMMDPARSER__HPP

#include "ctranscript.hpp"

BEGIN_OLIGOFAR_SCOPES

template<class TTranscript>
class CTranscriptWalker
{
public:
    typedef TTranscript transcript_type;
    typedef CTranscriptWalker<transcript_type> self_type;
    typedef typename TTranscript::const_iterator transcript_iterator;
    typedef typename TTranscript::data_type data_type;
    CTranscriptWalker( const transcript_type& t ) : m_transcript( t ), m_cursor( t.begin() ), m_step(0) {}
    self_type& operator ++ ();
    self_type  operator ++ ( int ) { self_type x( *this ); ++*this; return x; }
    self_type& operator += ( int );
    data_type  operator * () const;
    CTrBase::EEvent GetEvent() const { return (**this).GetEvent(); }
    int GetCount() const { return (**this).GetCount(); }
    bool       done() const;
    operator   bool () const { return !done(); }
protected:
    const transcript_type& m_transcript;
    transcript_iterator m_cursor;
    int m_step;
};

class CMdtagWalker
{
public:
    // seems like we don't care about deletion... it is known from CIGAR
    enum EEvent { eEvent_A = 'A', eEvent_C = 'C', eEvent_G = 'G', eEvent_T = 'T', eEvent_N = 'N', eEvent_equ = '=', eEvent_del = '^' /* optional */ };
    typedef CMdtagWalker self_type;
    CMdtagWalker( const char * tag ) : 
        m_cursor( tag ), m_step(0), 
        m_count( isdigit( *m_cursor ) ? strtol(m_cursor, const_cast<char**>(&m_cursor), 10) : ( ++m_cursor, 0L ) ) {
            if( m_count == 0 ) while( m_cursor[-1] == '^' ) ++m_cursor;
        } //: *tag == '^' ? ++m_cursor, 1 : *m_cursor == 0 ? 0 : 1 ) {}
    self_type& operator ++ ();
    self_type  operator ++ ( int ) { self_type x( *this ); ++*this; return x; }
    self_type& operator += ( int );
    int GetCount() const { return m_count ? m_count - m_step : 1; }
    int GetEvent() const { return m_count ? eEvent_equ : m_cursor[-1]; }
    operator bool () const { return GetEvent(); }
protected:
    const char * m_tag;
    const char * m_cursor;
    int  m_step;
    int  m_count;
};

template<class TDstCigar, class TSrcCigar>
inline TDstCigar& ParseSamMdTag( TDstCigar& dst, const TSrcCigar& src, const char * mdz ) 
{
    CTranscriptWalker<TSrcCigar> itrans( src );
    CMdtagWalker imdtag( mdz );

    while( imdtag && itrans ) {
        int step = min( imdtag.GetCount(), itrans.GetCount() );
        int mstep = step, tstep = step;
        switch( itrans.GetEvent() ) {
         case CTrBase::eEvent_SoftMask:
         case CTrBase::eEvent_HardMask:
         case CTrBase::eEvent_Padding:
         case CTrBase::eEvent_Insertion:
         case CTrBase::eEvent_Splice:
             mstep = 0; //dst.AppendItem( itrans.GetEvent(), step ); break;
         case CTrBase::eEvent_Overlap:
         case CTrBase::eEvent_Deletion:
             dst.AppendItem( itrans.GetEvent(), step );
             break;
         default:
             if( imdtag.GetEvent() == '=' ) {
                 dst.AppendItem( CTrBase::eEvent_Match, step );
             } else {
                 dst.AppendItem( CTrBase::eEvent_Replaced, step );
             }
             break;
        }
        imdtag += mstep;
        itrans += tstep;
    }

    if( imdtag ) {
        ostringstream err;
        while( imdtag ) {
            err << char(imdtag.GetEvent()) << "." << imdtag.GetCount() << ";";
            imdtag += imdtag.GetCount();
        }
        THROW( logic_error, "When parsing CIGAR of " << src.ToString() << " with MD:Z:" << mdz << " unprocessed MD:Z:" << err.str() << " left" );
    }
    while( itrans ) {
        switch( itrans.GetEvent() ) {
            case CTrBase::eEvent_SoftMask:
            case CTrBase::eEvent_HardMask:
            case CTrBase::eEvent_Insertion:
                dst.AppendItem( itrans.GetEvent(), itrans.GetCount() );
                itrans += itrans.GetCount();
                break;
            default:
                THROW( logic_error, "When parsing CIGAR of " << src.ToString() << " with MD:Z:" << mdz << " unprocessed CIGAR " << itrans.GetCount() << CTrBase::Event2Char( itrans.GetEvent() ) << " left" );
        }
    }
    return dst;
}

////////////////////////////////////////////////////////////////////////
// CTranscriptWalker implementation

template<class TTranscript>
inline CTranscriptWalker<TTranscript>& CTranscriptWalker<TTranscript>::operator ++ () 
{
    if( ++m_step >= m_cursor->GetCount() ) {
        ++m_cursor;
        m_step = 0;
    }
    return *this;
}

template<class TTranscript>
inline CTranscriptWalker<TTranscript>& CTranscriptWalker<TTranscript>::operator += ( int incr ) 
{
    while( m_step + incr >= m_cursor->GetCount() ) {
        incr -= (m_cursor->GetCount() - m_step);
        ++m_cursor;
        m_step = 0;
        if( m_cursor == m_transcript.end() ) { incr = 0; break; }
    }
    m_step += incr;
    return *this;
}

template<class TTranscript>
inline typename CTranscriptWalker<TTranscript>::data_type CTranscriptWalker<TTranscript>::operator * () const 
{
    return data_type( m_cursor->GetEvent(), m_cursor->GetCount() - m_step );
}

template<class TTranscript>
inline bool CTranscriptWalker<TTranscript>::done() const 
{
    return !(m_cursor != m_transcript.end());
}

////////////////////////////////////////////////////////////////////////
// CMdtagWalker implementation

inline CMdtagWalker& CMdtagWalker::operator ++ () 
{
    while( true ) {
        if( m_count == 0 ) {
            while( *m_cursor == '^' ) ++m_cursor;
            if( m_cursor[-1] == 0 ) return *this;
            else if( isdigit( *m_cursor ) ) {
                m_count = strtol( m_cursor, const_cast<char**>(&m_cursor), 10 );
                m_step = 0;
            } else m_cursor++;
            return *this;
        } else if( ++m_step >= m_count ) {
            m_step = 0;
            m_count = 0;
        } else return *this;
    }
}

inline CMdtagWalker& CMdtagWalker::operator += ( int i )
{
    while( i > 0 ) {
retry:
        if( m_count == 0 ) {
            while( *m_cursor == '^' ) ++m_cursor;
            if( m_cursor[-1] == 0 ) return *this;
            else if( isdigit( *m_cursor ) ) {
                m_count = strtol( m_cursor, const_cast<char**>(&m_cursor), 10 );
                m_step = 0;
                if( m_count > --i ) {
                    m_step  += i;
                    return *this;
                } else {
                    i -= m_count;
                    m_count = 0;
                    m_step = 0;
                }
            } else { if( !++m_cursor ) break; --i; }
        } else if( m_step + i >= m_count ) {
            i -= ( m_count - m_step );
            m_step = 0;
            m_count = 0;
            goto retry;
        } else { 
            m_step += i;
            return *this;
        }
    }
    return *this;
}


END_OLIGOFAR_SCOPES

#endif
