#ifndef OLIGOFAR_DUST__HPP
#define OLIGOFAR_DUST__HPP

#include "cringbuffer.hpp"
#include "cseqcoding.hpp"
#include <assert.h>

BEGIN_OLIGOFAR_SCOPES

class CComplexityMeasure
{
public:
    CComplexityMeasure( int windowSize = 2 );
    operator double() const { return Get(); }
    void Add( int tripletId );
    void Del( int tripletId );
    void IncWin() { ++m_windowSize; } // to take into account ambiguities in symmetrical way
    void DecWin() { --m_windowSize; }
    double Get() const {  return m_windowSize > 3 ? m_score / 2 / ( m_windowSize - 3 ) : 0; }
    int GetWindowSize() const { return m_windowSize; }
    void Reset();
    string Debug() const { ostringstream o; for( int i = 0; i < 64; ++i ) { if( i ) o << " "; o << setw(2) << m_tripletCounts[i]; } return o.str(); }
protected:
    double m_score;
    int m_windowSize;
    int m_tripletCounts[64];
};

class CComplexityMeasureWindow : public CComplexityMeasure
{
public: 
    CComplexityMeasureWindow() : CComplexityMeasure(0) {}

    void Del();
    void Add( const CIupacnaBase& b ) { Add( CNcbi8naBase( b ) ); }
    void Add( const CNcbi8naBase& b ) { Add( b.GetSmallestNcbi2na() ); }
    void Add( const CNcbi2naBase& );

    void Reset() { m_buffer.clear(); CComplexityMeasure::Reset(); }
private:
    void Add( int tripletId );
    void Del( int tripletId );
    void IncWin() {}
    void DecWin() {}
protected:
    CRingBuffer<char,true> m_buffer;
};

////////////////////////////////////////////////////////////////////////

inline CComplexityMeasure::CComplexityMeasure( int windowSize ) 
    : m_score( windowSize > 1 ? (windowSize - 2)*(windowSize - 3) : 0 ), m_windowSize( windowSize )
{
    fill( m_tripletCounts, m_tripletCounts + 64, 0 );
    if( m_windowSize > 2 ) m_tripletCounts[0] = m_windowSize - 2;
}

inline void CComplexityMeasure::Add( int tripletId ) 
{
    assert( tripletId >= 0 );
    assert( tripletId < 64 );
    
    if( m_windowSize >= 2 ) {
        int & tripletCount = m_tripletCounts[tripletId];
        m_score -= tripletCount * ( tripletCount - 1 );
        m_score += tripletCount * ( tripletCount + 1 );
        tripletCount ++;
    }
    m_windowSize ++;
}

inline void CComplexityMeasure::Del( int tripletId )
{
    if( m_windowSize > 2 ) {
        int & tripletCount = m_tripletCounts[tripletId];
        if( tripletCount <= 0 ) 
			throw logic_error( "CComplexityMeasure::Del(int tripletId): Deleting non-existing triplet!" );
        tripletCount --;
        m_score -= tripletCount * ( tripletCount + 1 );
        m_score += tripletCount * ( tripletCount - 1 );
    }
    m_windowSize -- ;
}

inline void CComplexityMeasure::Reset()
{
    m_windowSize = 2;
    m_score = 0;
    fill( m_tripletCounts, m_tripletCounts + 64, 0 );
}

///////////////////////////////////////////////////////////////////////
inline void CComplexityMeasureWindow::Add( const CNcbi2naBase& base ) 
{
    switch( GetWindowSize() ) {
        case 0: m_buffer.push_back( base ); break;
        case 1:
        case 2: (m_buffer.back() <<= 2) |= base; break;
        default: m_buffer.push_back( ((m_buffer.back() << 2) | base) & 0x3f ); break;
    }
    CComplexityMeasure::Add( m_buffer.back() );
}

inline void CComplexityMeasureWindow::Del() 
{ 
    int tripletId = 0;
    switch( GetWindowSize() ) {
        case 0: THROW( logic_error, "Deleting base from empty CComplexityMeasureWindow" );
        case 1: 
        case 2:
        case 3: tripletId = (m_buffer.front() >> (GetWindowSize()-1)); break;
        default: tripletId = m_buffer.pop_front(); break;
    }
    CComplexityMeasure::Del( tripletId ); 
}

END_OLIGOFAR_SCOPES

#endif
