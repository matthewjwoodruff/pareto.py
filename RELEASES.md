# Release Notes

## 1.1.1

Improvements to the API
* eps_sort can now handle all these kinds of input:
    - DataFrame
    - [DataFrame]
    - ndarray
    - [ndarray] 
    - generator 
    - [generator] 
    - generator of generators 
    - [double-subscriptable list of floats] 
    - double-subscriptable list of floats
    - [double-subscriptable list of strings]
    - double-subscriptable list of strings
* Added flag_nondominated function to return a list of booleans, or a 
  list of lists if it is called with a list of tables instead of a single 
  table

## 1.0 and before

Developed basic functionality.

