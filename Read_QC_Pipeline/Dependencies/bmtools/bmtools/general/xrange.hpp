#ifndef OLIGOFAR_XRANGE__HPP
#define OLIGOFAR_XRANGE__HPP

#include "of-debug.hpp"
#include <iostream>

BEGIN_OLIGOFAR_SCOPES

////////////////////////////////////////////////////////////////////////
// Written because CRange does something weird on SetFrom, SetTo
//
template<class valtype>
class CXRange
{
public:
    typedef valtype value_type;
    typedef CXRange<value_type> self_type;
    CXRange() : m_start(0), m_length(0) {}
    template<class T>
    CXRange( const T& other ) : m_start( other.GetFrom() ), m_length( other.GetLength() ) {}
    CXRange& SetOpen( value_type begin, value_type end ) { m_start = begin; m_length = end > begin ? end - begin : 0; return *this; }
    CXRange& SetClosed( value_type from, value_type to ) { m_start = from;  if( to > from ) m_length = to - from + value_type(1); else m_length = 0; return *this; }
    CXRange& SetOffLen( value_type from, value_type len ) { m_start = from; m_length = len; return *this; }
    static CXRange MakeOpen( value_type begin, value_type end ) { return CXRange().SetOpen( begin, end ); }
    static CXRange MakeClosed( value_type from, value_type to ) { return CXRange().SetClosed( from, to ); }
    static CXRange MakeOffLen( value_type from, value_type len ) { return CXRange().SetOffLen( from, len ); }
    self_type& SetFrom( value_type x ) { m_length -= (x - m_start ); m_start = x; return *this; }
    self_type& SetToOpen( value_type x ) { m_length = x - m_start; return *this; }
    self_type& SetTo( value_type x ) { m_length = x - m_start + 1; return *this; }
    self_type& SetLength( value_type x ) { m_length = x; return *this; }
    self_type& SetLengthDown( value_type x ) { m_start -= (x-m_length); m_length = x; return *this; }
    self_type& GrowBy( value_type x ) { m_start -= x; m_length += 2*x; return *this; }
    self_type& ShrinkBy( value_type x ) { m_start += x; m_length -= 2*x; return *this; }
    self_type& ShiftBy( value_type x ) { m_start += x; return *this; }
    self_type& ShiftTo( value_type x ) { m_start = x; return *this; }
    value_type GetLength() const { return m_length; }
    value_type GetFrom() const { return m_start; }
    value_type GetTo() const { return m_start + m_length - value_type(1); }
    value_type GetToOpen() const { return m_start + m_length; }
    value_type GetMiddle() const { return m_start + m_length/2; }
    value_type SetEmpty() { m_length = 0; }
    bool Empty() const { return m_length <= 0; }
    bool NotEmpty() const { return m_length > 0; }
    class incr { public: value_type operator ()( value_type& x ) const { return ++x; } };
    class decr { public: value_type operator ()( value_type& x ) const { return --x; } };
    template<class fwd, class rev>
    class iter
    {
    public:
        iter( value_type x ) : m_value(x) {}
        iter& operator ++ () { fwd m; m(m_value); return *this; }
        iter  operator ++ (int) { iter x(*this); ++*this; return x; }
        iter& operator -- () { rev m; m(m_value); return *this; }
        iter  operator -- (int) { iter x(*this); --*this; return x; }
        const value_type& operator *  () const { return  m_value; }
        const value_type* operator -> () const { return &m_value; }
        bool operator == ( const iter& other ) const { return m_value == other.m_value; }
        bool operator != ( const iter& other ) const { return m_value != other.m_value; }
    protected:
        value_type m_value;
    };
    typedef iter<incr,decr> const_iterator;
    typedef iter<decr,incr> const_reverse_iterator;
    const_iterator begin() const { return const_iterator( m_start ); }
    const_iterator end() const { return const_iterator( m_start + m_length ); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator( m_start + m_length - 1 ); }
    const_reverse_iterator rend() const { return const_reverse_iterator( m_start - 1 ); }
    template<class T>
    self_type& CombineWith( const CXRange<T>& other ) {
        return SetOpen( min( GetFrom(), other.GetFrom() ), max( GetToOpen(), other.GetToOpen() ) );
    }
    self_type CombinationWith( const self_type& other ) const {
        return MakeOpen( min( GetFrom(), other.GetFrom() ), max( GetToOpen(), other.GetToOpen() ) );
    }
    template<class T>
    self_type& IntersectWith( const CXRange<T>& other ) {
        return SetOpen( max( GetFrom(), other.GetFrom() ), min( GetToOpen(), other.GetToOpen() ) );
    }
    self_type IntersectionWith( const self_type& other ) const {
        return MakeOpen( max( GetFrom(), other.GetFrom() ), min( GetToOpen(), other.GetToOpen() ) );
    }
    template<class T>
    bool IntersectsWith( const CXRange<T>& other ) const {
        return IntersectingWith( other );
    }
    template<class T>
    bool IntersectingWith( const CXRange<T>& other ) const {
        return max( GetFrom(), other.GetFrom() ) < min( GetToOpen(), other.GetToOpen() );
    }
    template<class T>
    bool Contains( const CXRange<T>& other ) const {
        return GetFrom() <= other.GetFrom() && GetToOpen() >= other.GetToOpen();
    }
    template<class T>
    bool ContainedIn( const CXRange<T>& other ) const { return other.Contains( *this ); }
    template<class T>
    bool Contains( T pos ) const { return GetFrom() <= pos && GetTo() >= pos; }
    template<class T> self_type& operator &= ( const T& t ) { return IntersectWith( t ); }
    template<class T> self_type& operator |= ( const T& t ) { return CombineWith( t ); }
    template<class T> self_type& operator += ( const T& t ) { return CombineWith( t ); }
    /*
    friend self_type operator & ( const self_type& a, const self_type& b ) { return a.IntersectionWith(b); }
    friend self_type operator | ( const self_type& a, const self_type& b ) { return a.CombinationWith(b); }
    friend self_type operator + ( const self_type& a, const self_type& b ) { return a.CombinationWith(b); }
    */
    template<class T> bool operator == ( const T& other ) const { return GetFrom() == other.GetFrom() && GetLength() == other.GetLength(); }
    template<class T> bool operator != ( const T& other ) const { return GetFrom() != other.GetFrom() || GetLength() != other.GetLength(); }
    template<class T> bool operator <= ( const T& other ) const {
        if( GetFrom() < other.GetFrom() ) return true;
        if( GetFrom() > other.GetFrom() ) return false;
        return GetLength() <= other.GetLength();
    }
    template<class T> bool operator >= ( const T& other ) const {
        if( GetFrom() > other.GetFrom() ) return true;
        if( GetFrom() < other.GetFrom() ) return false;
        return GetLength() >= other.GetLength();
    }
    template<class T> bool operator <  ( const T& other ) const {
        if( GetFrom() < other.GetFrom() ) return true;
        if( GetFrom() > other.GetFrom() ) return false;
        return GetLength() <  other.GetLength();
    }
    template<class T> bool operator >  ( const T& other ) const {
        if( GetFrom() > other.GetFrom() ) return true;
        if( GetFrom() < other.GetFrom() ) return false;
        return GetLength() >  other.GetLength();
    }
protected:
    value_type m_start;
    value_type m_length;
};

END_OLIGOFAR_SCOPES

namespace std {
    template<class vtype>
    ostream& operator << ( ostream& o, const OLIGOFAR_SCOPES::CXRange<vtype> r ) {
        return o << "[" << r.GetFrom() << "+" << r.GetLength() << "=" << r.GetToOpen() << ")";
    }
}

#endif
