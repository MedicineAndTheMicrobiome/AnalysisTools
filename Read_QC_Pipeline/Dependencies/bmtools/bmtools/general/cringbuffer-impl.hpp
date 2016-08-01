#ifndef CRINGBUFFER_IMPL__HPP
#define CRINGBUFFER_IMPL__HPP

#ifndef OLIGOFAR__NOINCLUDES
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <assert.h>
#endif

template<class Item,bool simpleType = false>
class CRingBuffer
{
public:
    typedef size_t size_type;
    typedef Item value_type;
    typedef Item& ref_type;
    typedef Item* ptr_type;
    typedef const Item& const_ref_type;
    typedef const Item* const_ptr_type;
    typedef CRingBuffer<Item,simpleType> self_type;

    ~CRingBuffer();
    CRingBuffer( size_type reserve = 0 );

    explicit CRingBuffer( const self_type& other );
    self_type& operator = ( const self_type& other );

    bool empty() const { return m_size == 0; }
    size_type size() const { return m_size; }
    size_type avail() const { return m_capacity - m_size; }
    size_type capacity() const { return m_capacity; }

    void clear();
    void free();
    void reserve( size_type );
    ref_type back() { return m_data[(m_head + m_size - 1)%m_capacity]; }
    ref_type front() { return m_data[m_head]; }
    ref_type operator [] ( int i ) { return m_data[(m_head + i)%m_capacity]; }
    const_ref_type back() const { return m_data[(m_head + m_size - 1)%m_capacity]; }
    const_ref_type front() const { return m_data[m_head]; }
    const_ref_type operator [] ( int i ) const { return m_data[(m_head + i)%m_capacity]; }
    void push_back( const_ref_type item );
    void push_front( const_ref_type item );
    value_type pop_back();
    value_type pop_front();

private:
    static ptr_type x_copy( const_ptr_type begin, const_ptr_type end, ptr_type dest );
protected:
    value_type * m_data;
    size_type m_capacity;
    size_type m_size;
    size_type m_head;
};

////////////////////////////////////////////////////////////////////////
// Implementation

template<class Item,bool simpleType>
inline CRingBuffer<Item,simpleType>::~CRingBuffer() 
{
    delete[] m_data;
}

template<class Item,bool simpleType>
inline CRingBuffer<Item,simpleType>::CRingBuffer( size_type reserve ) :
    m_data(0),
    m_capacity(0),
    m_size(0),
    m_head(0)
{
    if( reserve ) m_data = new value_type[m_capacity = reserve];
}

template<class Item,bool simpleType>
inline CRingBuffer<Item,simpleType>::CRingBuffer( const self_type& other ) :
    m_data(new Item[other.m_capacity]),
    m_capacity( other.m_capacity ),
    m_size( other.m_size ),
    m_head( other.m_head )
{
    x_copy( other.m_data, other.m_data + m_capacity, m_data );
}

template<class Item,bool simpleType>
inline CRingBuffer<Item,simpleType>& CRingBuffer<Item,simpleType>::operator = ( const self_type& other )
{
    if( &other == this ) return *this;
    delete[] m_data;
    m_data = new Item[other.m_capacity];
    m_capacity = other.m_capacity;
    m_size = other.m_size;
    m_head = other.m_head;
    x_copy( other.m_data, other.m_data + m_capacity, m_data );
}

template<class Item,bool simpleType>
inline void CRingBuffer<Item,simpleType>::clear()
{
    m_head = m_size = 0;
}

template<class Item,bool simpleType>
inline void CRingBuffer<Item,simpleType>::free()
{
    clear();
    m_capacity = 0;
    delete[] m_data;
    m_data = 0;
}

template<class Item,bool simpleType>
inline Item * CRingBuffer<Item,simpleType>::x_copy( const_ptr_type begin, const_ptr_type end, ptr_type dest ) 
{
    assert( end >= begin );
    if( simpleType ) {
        size_type count = end - begin;
        if( count ) std::memmove( dest, begin, count*sizeof(value_type) );
        return dest + count;
    } else {
        return ::std::copy( begin, end, dest );
    }
}

template<class Item,bool simpleType>
inline void CRingBuffer<Item,simpleType>::reserve( size_type newsz ) 
{
    if( newsz < m_size ) return;
    value_type * new_data = new value_type[newsz];
    value_type * end = new_data;
    size_type tail = m_head + m_size;
    if( tail > m_capacity ) {
        end = x_copy( m_data + m_head, m_data + m_capacity       , end );
        end = x_copy( m_data         , m_data - m_capacity + tail, end );
    } else {
        end = x_copy( m_data + m_head, m_data + tail, end );
    }
    delete[] m_data;
    m_data = new_data;
    m_head = 0;
    m_capacity = newsz;
}

template<class Item,bool simpleType>
inline void CRingBuffer<Item,simpleType>::push_back( const_ref_type item )
{
    if( m_capacity == 0 ) reserve( 16 );
    if( !avail() ) reserve( capacity()*2 );
    m_data[(m_head + m_size)%m_capacity] = item;
    m_size++;
}

template<class Item,bool simpleType>
inline void CRingBuffer<Item,simpleType>::push_front( const_ref_type item )
{
    if( m_capacity == 0 ) reserve( 16 );
    if( !avail() ) reserve( capacity()*2 );
    if( m_head == 0 ) m_head = m_capacity;
    m_head--;
    m_size++;
    m_data[m_head] = item;
}

template<class Item,bool simpleType>
inline Item CRingBuffer<Item,simpleType>::pop_back()
{
    if( empty() ) throw std::range_error( "Failure to pop_back from empty CRingBuffer" );
    --m_size;
    return m_data[(m_head + m_size)%m_capacity];
}

template<class Item,bool simpleType>
inline Item CRingBuffer<Item,simpleType>::pop_front()
{
    if( empty() ) throw std::range_error( "Failure to pop_back from empty CRingBuffer" );
    size_type t = m_head;
    --m_size;
    ++m_head %= m_capacity;
    return m_data[t%m_capacity];
}

#endif
