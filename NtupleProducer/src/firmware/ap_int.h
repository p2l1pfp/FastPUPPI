template<typename T>
class bitref_ {
    public:
        bitref_(T & number, const int nbit) : num_(number), bit_(nbit) {}
        operator bool() const { return (num_ & (1 << bit_)); }
        bitref_ & operator=(bool value) {
            if (value) num_ |=  (1 << bit_);
            else       num_ &= ~(1 << bit_); 
            return *this;
        }
    private:
        T & num_;
        const int bit_;
};

template<unsigned int N, typename T=int> class ap_int {
    public:
        typedef T base_type;

        ap_int(T value=0) : value_(mask(value)) { /*printf("ap_int<%u>(%d) = { value_ = %d; }\n", N, value, value_);*/ }

        template<unsigned int M> 
        ap_int(const ap_int<M,T> & other) { value_ = mask(other.value_); }

        operator T() const { return value_; }

        static T mask(T value) { return (value >= 0 ? value & pos_mask() : value | neg_mask()); }
        static T pos_mask() { return (1 << (N-1))-1; }
        static T neg_mask() { return    ~pos_mask(); }

        ap_int<N,T> & operator+=(int i) { value_ = mask(value_ + i); return *this; }
        ap_int<N,T> & operator-=(int i) { value_ = mask(value_ - i); return *this; }

        bool operator[](int i) const { return value_ & (1 << i); }
        bitref_<T> operator[](int i) { return bitref_<T>(value_, i); }

        ap_int<N,T> operator<<(int i) const { return ap_int<N,T>(value_ << i); }
        ap_int<N,T> operator>>(int i) const { return ap_int<N,T>(value_ >> i); }
    private:
        T value_;
};

template<unsigned int N, unsigned int M, typename T>
inline int operator+ (const ap_int<N,T> & a, const ap_int<M,T> &b ) { return T(a)+T(b); }
template<unsigned int N, unsigned int M, typename T>
inline int operator- (const ap_int<N,T> & a, const ap_int<M,T> &b ) { return T(a)-T(b); }
template<unsigned int N, unsigned int M, typename T>
inline int operator* (const ap_int<N,T> & a, const ap_int<M,T> &b ) { return T(a)*T(b); }

template<unsigned int N, typename T=unsigned> class ap_uint {
    public:
        typedef T base_type;

        ap_uint(const T &value=0) : value_(value & mask()) {}

        template<unsigned int M>
        ap_uint(const ap_uint<M,T> & other) { value_ = (other & mask()); }

        operator T() const { return value_; }

        static T mask() { return (1 << N)-1; }

        ap_uint<N,T> & operator+=(unsigned i) { value_ = mask() & (value_ + i); return *this;  }
        ap_uint<N,T> & operator-=(unsigned i) { value_ = mask() & (value_ - i); return *this;  }

        bool operator[](int i) const { return value_ & (1 << i); }
        bitref_<T> operator[](int i) { return bitref_<T>(value_, i); }

        ap_uint<N,T> operator<<(int i) const { return ap_uint<N,T>(value_ << i); }
        ap_uint<N,T> operator>>(int i) const { return ap_uint<N,T>(value_ >> i); }
    private:
        T value_;
};

template<unsigned int N, unsigned int M, typename T>
inline int operator+ (const ap_uint<N,T> & a, const ap_uint<M,T> &b ) { return T(a)+T(b); }
template<unsigned int N, unsigned int M, typename T>
inline int operator- (const ap_uint<N,T> & a, const ap_uint<M,T> &b ) { return T(a)-T(b); }
template<unsigned int N, unsigned int M, typename T>
inline int operator* (const ap_uint<N,T> & a, const ap_uint<M,T> &b ) { return T(a)*T(b); }

