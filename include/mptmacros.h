/*****************************************************************************
 *
 *
 * @file   mptmacros.h
 * @author
 * @date
 *
 *****************************************************************************/
#ifndef MPTMACROS_H_
#  define MPTMACROS_H_ 1
/* ------------------------------------------------------------------------- */
const double D_PI = 3.14159265;
const double TWO_PI = D_PI * 2.0;
const double FOURTHIRD_PI = D_PI * 4.0 / 3.0;
const double FOUR_PI = D_PI * 4.0;

const double D_VOL_FRAC = D_PI*1.0E-6;
const double D_D_CONST = 3E-9;

const double D_GI = 2.6752218744E8;//(radsT)^-1 //NIST (gyromagnetic ratio H^1)


//calculation of constant part of magnetic field
const double D_MO = FOUR_PI * 1.0E-7; // μ0 (Hardy)

const double D_M_IONP = 3.8E5; //Magnetization

const double B_fixed = D_MO*D_M_IONP*(1./3.); //equatorial field as defined in Vuong paper

const double D_T_STEP = 10E-6;// duration of step (10 µs original)

template <class T> //instead of overload
T square(T a) {
    return a * a;
}

template <class T>
T cube(T a) {
    return a * a * a;
}
/* ------------------------------------------------------------------------- */
/** @name CAST macros                                                        */
/**@{*/
/** Macro to static cast a value to double type                              */
#  define SC_D( VAL ) ( static_cast< double >( VAL ) )
/** Macro to static cast a value to a long type                              */
#  define SC_L( VAL ) ( static_cast< long >( VAL ) )
/** Macro to static cast a value to a int32_t type                           */
#  define SC_I32( VAL ) ( static_cast< int32_t >( VAL ) )
/** Macro to static cast a value to a int16_t type                           */
#  define SC_I16( VAL ) ( static_cast< int16_t >( VAL ) )
/** Macro to static cast a value to a unsigned short type                    */
#  define SC_US( VAL ) ( static_cast< unsigned short >( VAL ) )
/** Macro to static cast a value to a unsigned int type                      */
#  define SC_UI( VAL ) ( static_cast< unsigned int >( VAL ) )
/** Macro to static cast a value to a int type                               */
#  define SC_I( VAL ) ( static_cast< int >( VAL ) )
/** Macro to static cast a value to a unsigned long (int) type               */
#  define SC_UL( VAL ) ( static_cast< unsigned long >( VAL ) )
/** Macro to static cast a value to a std::string type                       */
#  define SC_S( VAL ) ( static_cast< std::string >( VAL ) )
/**@}*/

/* ------------------------------------------------------------------------- */
#endif // #ifndef MPTMACROS_H_
/* ---- EOF ---------------------------------------------------------------- */
