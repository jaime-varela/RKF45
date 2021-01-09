#include <concepts>


namespace RungeKutta
{
  template<class T>
  concept NumT = std::is_arithmetic<T>::value; 

  template<class T>
  concept Index = std::is_integral<T>::value;

  template<class container,class number>
  concept CoordinateContainer = std::is_arithmetic<number>::value &&
  requires(container A,number b) {
    {A[0] < A[0]} -> std::same_as<bool>;
    {A[0] + b} -> std::same_as<number>;
  };

}