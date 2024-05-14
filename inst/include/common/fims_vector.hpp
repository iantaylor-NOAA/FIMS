#ifndef FIMS_VECTOR_HPP
#define FIMS_VECTOR_HPP

#include "../interface/interface.hpp"

#define FIMS_VECTOR_IMLICIT_THROW
#define FIMS_FULL_VECTOR_OPERATORS

namespace fims {

    /**
     * Wrapper class for std::vector types. If this file is compiled with
     * -DTMB_MODEL, conversion operators are defined for TMB vector types.
     *
     * All std::vector functions are copied over from the std library. While some of
     * these may not be called explicitly in FIMS, they may be required to run other
     * std library functions.
     *
     */
    template <typename Type>
    class Vector {
        std::vector<Type> vec_m;

        /**
         * @brief friend comparison operator. Allows the operartor to see private
         * members of fims::Vector<Type>.
         */
        template <typename T>
        friend bool operator==(const fims::Vector<T>& lhs,
                const fims::Vector<T>& rhs);

    public:
        // Member Types

        typedef
        typename std::vector<Type>::value_type value_type; /*!<Member type Type>*/
        typedef typename std::vector<Type>::allocator_type
        allocator_type; /*!<Allocator for type Type>*/
        typedef typename std::vector<Type>::size_type size_type; /*!<Size type>*/
        typedef typename std::vector<Type>::difference_type
        difference_type; /*!<Difference type>*/
        typedef typename std::vector<Type>::reference
        reference; /*!<Reference type &Type>*/
        typedef typename std::vector<Type>::const_reference
        const_reference; /*!<Constant eference type const &Type>*/
        typedef typename std::vector<Type>::pointer pointer; /*!<Pointer type Type*>*/
        typedef typename std::vector<Type>::const_pointer
        const_pointer; /*!<Constant ointer type const Type*>*/
        typedef typename std::vector<Type>::iterator iterator; /*!<Iterator>*/
        typedef typename std::vector<Type>::const_iterator
        const_iterator; /*!<Constant iterator>*/
        typedef typename std::vector<Type>::reverse_iterator
        reverse_iterator; /*!<Reverse iterator>*/
        typedef typename std::vector<Type>::const_reverse_iterator
        const_reverse_iterator; /*!<Constant reverse iterator>*/

#ifdef FIMS_FULL_VECTOR_OPERATORS

        //friend operators
        template<typename T>
        friend Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2);
        template<typename T>
        friend Vector<T> operator+(const Vector<T>& v1, const T& a);
        template<typename T>
        friend Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2);
        template<typename T>
        friend Vector<T> operator-(const Vector<T>& v1, const T& a);
        template<typename T>
        friend Vector<T> operator*(const Vector<T>& v1, const Vector<T>& v2);
        template<typename T>
        friend Vector<T> operator*(const Vector<T>& v, const T& a);
        template<typename T>
        friend Vector<T> operator*(const T& a, const Vector<T>& v);
        template<typename T>
        friend Vector<T> operator/(const Vector<T>& v, const T& a); 
#endif
        // Constructors

        /**
         * Default constructor.
         */
        Vector() {
            this->vec_m.resize(1);
        }

        Vector(const std::initializer_list<Type>& l) {
            this->vec_m = std::vector<Type>(l);
        }

        /**
         * @brief Constructs a Vector of length "size" and sets the elements with the
         * value from input "value".
         */
        Vector(size_t size, const Type& value = Type()) {
            this->vec_m.resize(size, value);
        }

        /**
         * @brief Copy constructor.
         */
        Vector(const Vector<Type>& other) {
            this->vec_m.resize(other.size());
            for (size_t i = 0; i < this->vec_m.size(); i++) {
                this->vec_m[i] = other[i];
            }
        }

        /**
         * @brief Initialization constructor from std::vector<Type> type.
         */
        Vector(const std::vector<Type>& other) {
            this->vec_m = other;
        }

        // TMB specific constructor
#ifdef TMB_MODEL

        /**
         * @brief Initialization constructor from tmbutils::vector<Type> type.
         */
        Vector(const tmbutils::vector<Type>& other) {
            this->vec_m.resize(other.size());
            for (size_t i = 0; i < this->vec_m.size(); i++) {
                this->vec_m[i] = other[i];
            }
        }

#endif


        /**
         * The following are std::vector functions copied over from the standard
         * library. While some of these may not be called explicitly in FIMS, they may
         * be required to run other std library functions.
         */

        /**
         * @brief Returns a reference to the element at specified location pos. No
         * bounds checking is performed.
         */
        inline Type& operator[](size_t pos) {
            return this->vec_m[pos];
        }

        /**
         * @brief Returns a constant  reference to the element at specified location
         * pos. No bounds checking is performed.
         */
        inline const Type& operator[](size_t n) const {
            return this->vec_m[n];
        }

        /**
         * @brief Returns a reference to the element at specified location pos. Bounds
         * checking is performed.
         */
        inline Type& at(size_t n) {
            return this->vec_m.at(n);
        }

        /**
         * @brief Returns a constant reference to the element at specified location
         * pos. Bounds checking is performed.
         */
        inline const Type& at(size_t n) const {
            return this->vec_m.at(n);
        }

        /**
         * @brief  Returns a reference to the first element in the container.
         */
        inline reference front() {
            return this->vec_m.front();
        }

        /**
         * @brief  Returns a constant reference to the first element in the container.
         */
        inline const_reference front() const {
            return this->vec_m.front();
        }

        /**
         * @brief  Returns a reference to the last element in the container.
         */
        inline reference back() {
            return this->vec_m.back();
        }

        /**
         * @brief  Returns a constant reference to the last element in the container.
         */
        inline const_reference back() const {
            return this->vec_m.back();
        }

        /**
         * @brief Returns a pointer to the underlying data array.
         */
        inline pointer data() {
            return this->vec_m.data();
        }

        /**
         * @brief Returns a constant pointer to the underlying data array.
         */
        inline const_pointer data() const {
            return this->vec_m.data();
        }

        // iterators

        /**
         * @brief Returns an iterator to the first element of the vector.
         */
        inline iterator begin() {
            return this->vec_m.begin();
        }

        /**
         * @brief Returns an iterator to the element following the last element of the
         * vector.
         */
        inline iterator end() {
            return this->vec_m.end();
        }

        /**
         * @brief Returns a reverse iterator to the first element of the reversed
         * vector. It corresponds to the last element of the non-reversed vector.
         */
        inline reverse_iterator rbegin() {
            return this->vec_m.rbegin();
        }

        /**
         * @brief Returns a reverse iterator to the element following the last element
         * of the reversed vector. It corresponds to the element preceding the first
         * element of the non-reversed vector.
         */
        inline reverse_iterator rend() {
            return this->vec_m.rend();
        }

        /**
         * @brief Returns a constant reverse iterator to the first element of the
         * reversed vector. It corresponds to the last element of the non-reversed
         * vector.
         */
        inline const_reverse_iterator rbegin() const {
            return this->vec_m.rbegin();
        }

        /**
         * @brief Returns a constant reverse iterator to the element following the
         * last element of the reversed vector. It corresponds to the element
         * preceding the first element of the non-reversed vector.
         */
        inline const_reverse_iterator rend() const {
            return this->vec_m.rend();
        }

        // capacity

        /**
         * @brief Checks whether the container is empty.
         */
        inline bool empty() {
            return this->vec_m.empty();
        }

        /**
         * @brief Returns the number of elements.
         */
        inline size_type size() const {
            return this->vec_m.size();
        }

        /**
         * @brief Returns the maximum possible number of elements.
         */
        inline size_type max_size() const {
            return this->vec_m.max_size();
        }

        /**
         * @brief Reserves storage.
         */
        inline void reserve(size_type cap) {
            this->vec_m.reserve(cap);
        }

        /**
         * @brief Returns the number of elements that can be held in currently
         * allocated storage.
         */
        inline size_type capacity() {
            return this->vec_m.capacity();
        }

        /**
         *  @brief Reduces memory usage by freeing unused memory.
         */
        inline void shrink_to_fit() {
            this->vec_m.shrink_to_fit();
        }

        // modifiers

        /**
         * @brief Clears the contents.
         */
        inline void clear() {
            this->vec_m.clear();
        }

        /**
         * @brief Inserts value before pos.
         */
        inline iterator insert(const_iterator pos, const Type& value) {
            return this->vec_m.insert(pos, value);
        }

        /**
         * @brief Inserts count copies of the value before pos.
         */
        inline iterator insert(const_iterator pos, size_type count,
                const Type& value) {
            return this->vec_m.insert(pos, count, value);
        }

        /**
         * @brief Inserts elements from range [first, last) before pos.
         */
        template <class InputIt>
        iterator insert(const_iterator pos, InputIt first, InputIt last) {
            return this->vec_m.insert(pos, first, last);
        }

        /**
         * @brief Inserts elements from initializer list ilist before pos.
         */

        iterator insert(const_iterator pos, std::initializer_list<Type> ilist) {
            return this->vec_m.insert(pos, ilist);
        }

        /**
         * @brief Constructs element in-place.
         */
        template <class... Args>
        iterator emplace(const_iterator pos, Args&&... args) {
            return this->vec_m.emplace(pos, std::forward<Args>(args)...);
        }

        /**
         * @brief Removes the element at pos.
         */
        inline iterator erase(iterator pos) {
            return this->vec_m.erase(pos);
        }

        /**
         * @brief Removes the elements in the range [first, last).
         */
        inline iterator erase(iterator first, iterator last) {
            return this->vec_m.erase(first, last);
        }

        /**
         * @brief Adds an element to the end.
         */
        inline void push_back(const Type&& value) {
            this->vec_m.push_back(value);
        }

        /**
         * @brief Constructs an element in-place at the end.
         */
        template <class... Args>
        void emplace_back(Args&&... args) {
            this->vec_m.emplace_back(std::forward<Args>(args)...);
        }

        /**
         * @brief Removes the last element.
         */
        inline void pop_back() {
            this->vec_m.pop_back();
        }

        /**
         * @brief Changes the number of elements stored.
         */
        inline void resize(size_t s) {
            this->vec_m.resize(s);
        }

        /**
         * @brief Swaps the contents.
         */
        inline void swap(Vector& other) {
            this->vec_m.swap(other.vec_m);
        }

        // end std::vector functions

        /**
         * Conversion operators
         */

        /**
         * @brief Converts fims::Vector<Type> to std::vector<Type>
         */
        inline operator std::vector<Type>() {
            return this->vec_m;
        }


        operator Type() {
#ifdef FIMS_VECTOR_IMLICIT_THROW
            if (this->vec_m.size()  == 0 || this->vec_m.size() > 1) {
                throw std::range_error("Implicit conversion from fims::Vector to type \"Type\", fims::Vector has size not equal to 1. No scalar representation available");
            }
#endif
            return this->vec_m[0];
        }

#ifdef TMB_MODEL

        /**
         * @brief Converts fims::Vector<Type> to tmbutils::vector<Type>const
         */
        operator tmbutils::vector<Type>() const {
            tmbutils::vector<Type> ret;
            ret.resize(this->vec_m.size());
            for (size_t i = 0; i < this->vec_m.size(); i++) {
                ret[i] = this->vec_m[i];
            }
            return ret;
        }

        /**
         * @brief Converts fims::Vector<Type> to tmbutils::vector<Type>
         */
        operator tmbutils::vector<Type>() {
            tmbutils::vector<Type> ret;
            ret.resize(this->vec_m.size());
            for (size_t i = 0; i < this->vec_m.size(); i++) {
                ret[i] = this->vec_m[i];
            }
            return ret;
        }

#endif

        Vector<Type>& operator++() {
            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] = this->vec_m[i] + static_cast<Type> (1.0);
            }
            return *this;
        }

        Vector<Type>& operator--() {
            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] = this->vec_m[i] - static_cast<Type> (1.0);
            }
            return *this;
        }

        Vector<Type>& operator++(int) {
            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] = this->vec_m[i] + static_cast<Type> (1.0);
            }
            return *this;
        }

        Vector<Type>& operator--(int) {
            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] = this->vec_m[i] - static_cast<Type> (1.0);
            }
            return *this;
        }
#ifdef FIMS_FULL_VECTOR_OPERATORS
        Vector<Type>& operator+=(const Vector<Type>& v) {

            if (this->vec_m.size() < v.size()) {
                this->resize(v.size());
            }

            for (size_t i = 0; i < v.size(); i++) {
                this->vec_m[i] += v[i];
            }
            return *this;
        }

        Vector<Type>& operator+=(const Type& a) {

            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] += a;
            }
            return *this;
        }

        Vector<Type>& operator-=(const Vector<Type>& v) {

            if (this->vec_m.size() < v.size()) {
                this->resize(v.size());
            }

            for (size_t i = 0; i < v.size(); i++) {
                this->vec_m[i] -= v[i];
            }
            return *this;
        }

        Vector<Type>& operator-=(const Type& a) {

            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] -= a;
            }
            return *this;
        }

        Vector<Type>& operator*=(const Vector<Type>& v) {

            if (this->vec_m.size() < v.size()) {
                this->resize(v.size());
            }

            for (size_t i = 0; i < v.size(); i++) {
                this->vec_m[i] *= v[i];
            }
            return *this;
        }

        Vector<Type>& operator*=(const Type& a) {

            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] *= a;
            }
            return *this;
        }

        Vector<Type>& operator/=(const Vector<Type>& v) {

            if (this->vec_m.size() < v.size()) {
                this->resize(v.size());
            }

            for (size_t i = 0; i < v.size(); i++) {
                this->vec_m[i] /= v[i];
            }
            return *this;
        }

        Vector<Type>& operator/=(const Type& a) {

            for (size_t i = 0; i < this->size(); i++) {
                this->vec_m[i] /= a;
            }
            return *this;
        }
#endif

    private:
    }; // end fims::Vector class

    /**
     * @brief Comparison operator.
     */
    template <class T>
    bool operator==(const fims::Vector<T>& lhs, const fims::Vector<T>& rhs) {
        return lhs.vec_m == rhs.vec_m;
    }




#ifdef FIMS_FULL_VECTOR_OPERATORS
    //binary operators

    template<typename T>
    fims::Vector<T> operator+(const fims::Vector<T>& v1, const fims::Vector<T>& v2) {
        size_t n;

        //  set n to be the length of the longest vector and create a vector
        //  of that length to be returned
        if (v1.size() > v2.size()) {
            n = v1.size();
        } else {
            n = v2.size();
        }

        fims::Vector<T> w(n);

        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] + v2.vec_m[i];
            }
        } else if (v1.size() > v2.size()) {
            for (int i = 0; i < v2.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] + v2.vec_m[i];
            }
            for (int i = v2.size(); i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i];
            }
            FIMS_LOG << "vector add - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";
        } else {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] + v2.vec_m[i];
            }
            for (int i = v1.size(); i < v2.size(); i++) {
                w.vec_m[i] = v2.vec_m[i];
            }
            FIMS_LOG << "vector add - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";
        }
        return w;


    }

    template<typename T>
    fims::Vector<T> operator+(const fims::Vector<T>& v, const T& a) {

        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = v.vec_m[i] + a;
        }

        return w;


    }

    template<typename T>
    fims::Vector<T> operator+(const T& a, const fims::Vector<T>& v) {

        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = a + v.vec_m[i];
        }

        return w;


    }

    template<typename T>
    fims::Vector<T> operator-(const fims::Vector<T>& v1, const fims::Vector<T>& v2) {
        size_t n;

        //  set n to be the length of the longest vector and create a vector
        //  of that length to be returned
        if (v1.size() > v2.size()) {
            n = v1.size();
        } else {
            n = v2.size();
        }

        fims::Vector<T> w(n);

        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] - v2.vec_m[i];
            }
        } else if (v1.size() > v2.size()) {
            for (int i = 0; i < v2.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] - v2.vec_m[i];
            }
            for (int i = v2.size(); i < v1.size(); i++) {
                w.vec_m[i] = -1.0 * v1.vec_m[i];
            }
            FIMS_LOG << "vector subtract - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";
        } else {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] - v2.vec_m[i];
            }
            for (int i = v1.size(); i < v2.size(); i++) {
                w.vec_m[i] = -1.0 * v2.vec_m[i];
            }
            FIMS_LOG << "vector subtract - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";
        }
        return w;


    }

    template<typename T>
    fims::Vector<T> operator-(const fims::Vector<T>& v, const T& a) {
        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = v.vec_m[i] - a;
        }

        return w;
    }

    template<typename T>
    fims::Vector<T> operator-(const T& a, const fims::Vector<T>& v) {
        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = a - v.vec_m[i];
        }

        return w;
    }

    template<typename T>
    fims::Vector<T> operator*(const fims::Vector<T>& v1, const fims::Vector<T>& v2) {
        size_t n;

        //  set n to be the length of the longest vector and create a vector
        //  of that length to be returned
        if (v1.size() > v2.size()) {
            n = v1.size();
        } else {
            n = v2.size();
        }

        fims::Vector<T> w(n, T());

        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] * v2.vec_m[i];
            }
        } else if (v1.size() > v2.size()) {
            for (int i = 0; i < v2.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] * v2.vec_m[i];
            }
            FIMS_LOG << "vector multiply - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";
        } else {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] * v2.vec_m[i];
            }
            FIMS_LOG << "vector multiply - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";
        }

        return w;

    }

    template<typename T >
    fims::Vector<T> operator*(const fims::Vector<T>& v, const T & a) {
        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = v.vec_m[i] * a;
        }

        return w;
    }

    template<typename T >
    fims::Vector<T> operator*(const T& a, const fims::Vector<T>& v) {
        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = v.vec_m[i] * a;
        }

        return w;
    }


template<typename T>
    fims::Vector<T> operator/(const fims::Vector<T>& v1, const fims::Vector<T>& v2) {
        size_t n;

        //  set n to be the length of the longest vector and create a vector
        //  of that length to be returned
        if (v1.size() > v2.size()) {
            n = v1.size();
        } else {
            n = v2.size();
        }

        fims::Vector<T> w(n, T());

        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] / v2.vec_m[i];
            }
        } else if (v1.size() > v2.size()) {
            for (int i = 0; i < v2.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] / v2.vec_m[i];
            }

            FIMS_LOG << "vector divide - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";

        } else {
            for (int i = 0; i < v1.size(); i++) {
                w.vec_m[i] = v1.vec_m[i] / v2.vec_m[i];
            }

            FIMS_LOG << "vector divide - vectors different lengths\n";
            FIMS_LOG << "extra entries of shorter vector assumed to be 0.0\n";
        }

        return w;

    }

    template<typename T >
    fims::Vector<T> operator/(const fims::Vector<T>& v, const T & a) {
        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = v.vec_m[i] / a;
        }

        return w;
    }

    template<typename T>
    fims::Vector<T> operator/(const T& a, const fims::Vector<T>& v) {
        fims::Vector<T> w(v.size());

        for (int i = 0; i < v.size(); i++) {
            w.vec_m[i] = a / v.vec_m[i];
        }

        return w;
    }
   
#endif

} // namespace fims

#endif
