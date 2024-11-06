#ifndef H_DATA_ACCESSING
#define H_DATA_ACCESSING

// Template function to access 3D array elements
template<typename T>
inline T _array_at_dim3(T* array, int i0, int i1, int i2, int dim0, int dim1, int dim2) {
    return *(array + i0 * dim1 * dim2 + i1 * dim2 + i2);
}

// Template function to access 2D array elements
template<typename T>
inline T _array_at_dim3(T* array, int i0, int i1, int dim0, int dim1) {
    return *(array + i0 * dim1 + i1);
}

// Template function to set value in 3D array
template<typename T>
inline void _array_set_dim3(T* array, T value, int i0, int i1, int i2, int dim0, int dim1, int dim2) {
    *(array + i0 * dim1 * dim2 + i1 * dim2 + i2) = value;
}

// Template function to set value in 2D array
template<typename T>
inline void _array_set_dim3(T* array, T value, int i0, int i1, int dim0, int dim1) {
    *(array + i0 * dim1 + i1) = value;
}

inline double alpha_get(double* alpha, int sym_id, int span_from, int span_to, int MS){
    return _array_at_dim3(alpha, sym_id, span_from, span_to, 0, MS, MS);
}

inline double beta_get(double* beta, int sym_id, int span_from, int span_to, int MS){
    return _array_at_dim3(beta, sym_id, span_from, span_to, 0, MS, MS);
}

inline void alpha_set(double* alpha, double value, int sym_id, int span_from, int span_to, int MS){
    _array_set_dim3(alpha, value, sym_id, span_from, span_to, 0, MS, MS);
}

inline void beta_set(double* beta, double value, int sym_id, int span_from, int span_to, int MS){
    _array_set_dim3(beta, value, sym_id, span_from, span_to, 0, MS, MS);
}

inline void alpha_increase(double* alpha, double increase_value, int sym_id, int span_from, int span_to, int MS){
    alpha[sym_id * MS * MS + span_from * MS + span_to] += increase_value;
}

inline void beta_increase(double* beta, double increase_value, int sym_id, int span_from, int span_to, int MS){
    beta[sym_id * MS * MS + span_from * MS + span_to] += increase_value;
}

#endif
