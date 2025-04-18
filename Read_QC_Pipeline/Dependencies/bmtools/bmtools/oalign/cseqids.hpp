#ifndef OLIGOFAR_CSEQIDS__HPP
#define OLIGOFAR_CSEQIDS__HPP

#include "of-scopes.hpp"
#ifndef STANDALONE
#include <objects/seqloc/Seq_id.hpp>
#endif
#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>

BEGIN_OLIGOFAR_SCOPES

class CSeqIds 
{
public:
#ifndef STANDALONE
	typedef objects::CBioseq::TId TIds;
	typedef objects::CSeq_id TSeqId;
#endif
    typedef list<string> TIds;
    typedef string TSeqId;

	class CSeqDef
	{
	public:
		typedef set<string> TIdSet;
		CSeqDef( int oid = -1 ) : m_oid( oid ) {}
		CSeqDef( int oid, const TIds& ids ) : m_oid( oid ) { AddIds( ids ); }
		const TIdSet& GetIdStrings() const { return m_ids; }
		string GetBestIdString() const;
		string GetFullIdString() const;
        Int8 GetGi() const;
#ifndef STANDALONE
		bool HasId( const TSeqId& seqId ) const { return HasId( seqId.AsFastaString() ); }
		void AddId( const TSeqId& seqId ) { AddId( seqId.AsFastaString() ); }
		void AddIds( const string& seqIds );
		bool HasAnyId( const string& ids ) const;
#endif
		bool HasId( const string&  seqId ) const { return m_ids.find( seqId ) != m_ids.end(); }
		void AddId( const string&  seqId ) { m_ids.insert( seqId ); }
		void AddIds( const TIds& ids );
		bool HasAnyId( const TIds& ids ) const;
		bool HasOid() const { return m_oid >= 0; }
		bool SetOid( int oid ); // returns false if sequence is already seen
		int  GetOid() const;
        int  GetLength() const { return m_length; }
        void SetLength( int l ) { m_length = l; }
	protected:
		bool m_oidSet;
		int m_oid;
        int m_length;
		TIdSet m_ids;
	};

	typedef map<string,int> TSeqId2Ord;
	typedef deque<CSeqDef> TOrd2SeqIds;
	typedef map<int,int> TOid2Ord;

	CSeqIds() : m_lastOrd( -1 ) {}

	void Clear() { m_oid2ord.clear(); m_seqid2ord.clear(); m_ord2seqids.clear(); m_lastOrd = -1; }
	bool IsEmpty() const { return m_ord2seqids.empty(); }
	int  Register( const TIds& ids, int oid = -1 );
	int  Register( const string& ids, int oid = -1 );
	int  Register( int oid ) { TIds ids; return Register( ids, oid ); }
	const CSeqDef& GetSeqDef( int ord ) const;
    void SetSeqLength( int ord, int length ) { m_ord2seqids[ord].SetLength( length ); }

	int GetOrdById( const TIds& ids ) const;
	int GetOrdByOid( int oid ) const;

	// following function is used to avoid calling GetOrdById after Register in independend modules
	int GetLastOrd() const { ASSERT( m_lastOrd != -1 ); return m_lastOrd; }
protected:
	TOrd2SeqIds m_ord2seqids;
	TSeqId2Ord  m_seqid2ord;
	TOid2Ord    m_oid2ord;

	int m_lastOrd;
};

inline Int8 CSeqIds::CSeqDef::GetGi() const
{
    ITERATE( TIdSet, id, m_ids ) {
        if( id->substr( 0, 3 ) == "gi|" ) 
            return FromStr<Int8>( id->substr( 3 ) );
//            return NStr::StringToInt8( id->substr( 3 ) );
    }
    return -1;
}
        

inline string CSeqIds::CSeqDef::GetBestIdString() const
{
	if( m_ids.size() ) return *m_ids.begin();
	else return "UNKNOWN";
}

inline string CSeqIds::CSeqDef::GetFullIdString() const
{
	if( m_ids.size() ) {
		string out;
		ITERATE( TIdSet, id, m_ids ) {
			if( id != m_ids.begin() ) out += "|";
			out += *id;
		}
        return out;
	} else return "UNKNOWN";
}

#ifndef STANDALONE
inline void CSeqIds::CSeqDef::AddIds( const string& seqIds )
{
	TIds ids;
	objects::CSeq_id::ParseFastaIds( ids, seqIds );
	AddIds( ids );
}
#endif

inline void CSeqIds::CSeqDef::AddIds( const TIds& ids )
{
#ifdef STANDALONE
	ITERATE( TIds, id, ids ) AddId( *id );
#else
	ITERATE( TIds, id, ids ) AddId( **id );
#endif
}

inline bool CSeqIds::CSeqDef::HasAnyId( const TIds& ids ) const 
{
#ifdef STANDALONE
	ITERATE( TIds, id, ids ) if( HasId( *id ) ) return true;
#else
	ITERATE( TIds, id, ids ) if( HasId( **id ) ) return true;
#endif
	return false;
}

#ifndef STANDALONE
inline bool CSeqIds::CSeqDef::HasAnyId( const string& seqIds ) const 
{
	TIds ids;
	objects::CSeq_id::ParseFastaIds( ids, seqIds );
	return HasAnyId( ids );
}
#endif

inline bool CSeqIds::CSeqDef::SetOid( int oid )
{
	if( oid < 0 )
		THROW( logic_error, "Attempt to set invalid oid " << oid << " for " << GetFullIdString() );
	if( HasOid() && oid != m_oid )
		//THROW( logic_error, "Can't redefine oid for " << GetFullIdString() << " from " << m_oid << " to " << oid );
        cerr << "WARNING: reference sequence id " << GetFullIdString() << " appears multiple times (oids: " << m_oid << ", " << oid << ")\n";
        return false;
	m_oid = oid; 
    return true;
}

inline int CSeqIds::CSeqDef::GetOid() const
{
	if( !HasOid() ) 
		THROW( logic_error, "Attempt to get oid when it is not set for " << GetFullIdString() );
	return m_oid;
}

inline int CSeqIds::Register( const TIds& ids, int oid )
{
	ITERATE( TIds, i, ids ) {
#ifdef STANDALONE
        const string& id = *i;
#else
		string id( (*i)->AsFastaString() );
#endif
		TSeqId2Ord::iterator x = m_seqid2ord.find( id );
		if( x != m_seqid2ord.end() ) {
			int ord = x->second;
			m_ord2seqids[ord].AddIds( ids );
			if( oid >= 0 ) {
				m_ord2seqids[ord].SetOid( oid );
				m_oid2ord.insert( make_pair( oid, ord ) );
			}
			return m_lastOrd = ord;
		}
	}
	if( oid >= 0 ) {
		TOid2Ord::iterator i = m_oid2ord.find( oid );
		if( i != m_oid2ord.end() ) {
			int ord = i->second;
			m_ord2seqids[ord].AddIds( ids );
			return m_lastOrd = ord;
		}
	}
	ASSERT( ids.size() || oid >= 0 );
	int ord = m_ord2seqids.size();
	m_ord2seqids.push_back( CSeqDef( oid, ids ) );
#ifdef STANDALONE
	ITERATE( TIds, i, ids ) 
		m_seqid2ord.insert( make_pair( *i, ord ) );
#else
	ITERATE( TIds, i, ids ) 
		m_seqid2ord.insert( make_pair( (*i)->AsFastaString(), ord ) );
#endif
	if( oid >= 0 ) m_oid2ord.insert( make_pair( oid, ord ) );
	return m_lastOrd = ord;
}

inline int CSeqIds::Register( const string& idstr, int oid )
{
	TIds ids;
#ifdef STANDALONE
    ids.push_back( idstr );
#else
	objects::CSeq_id::ParseFastaIds( ids, idstr );
#endif
	return Register( ids, oid );
}
	
inline const CSeqIds::CSeqDef& CSeqIds::GetSeqDef( int ord ) const
{
	ASSERT( ord >= 0 && ord < (int)m_ord2seqids.size() );
	return m_ord2seqids[ord];
}

inline int CSeqIds::GetOrdById( const TIds& ids ) const
{
	ITERATE( TIds, i, ids ) {
#ifdef STANDALONE
        const string& id = *i;
#else
		string id( (*i)->AsFastaString() );
#endif
		TSeqId2Ord::const_iterator x = m_seqid2ord.find( id );
		if( x != m_seqid2ord.end() ) return x->second;
	}
	return -1;
}

inline int CSeqIds::GetOrdByOid( int oid ) const
{
	TOid2Ord::const_iterator x = m_oid2ord.find( oid );
	if( x == m_oid2ord.end() ) return -1;
	return x->second;
}

END_OLIGOFAR_SCOPES

#endif
