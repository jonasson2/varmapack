// Define macros for convenient debug printing (see also debugprint-demo.c).
#ifndef DEBUGPRINT_H
#define DEBUGPRINT_H

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define printf mexPrintf
#endif

#ifdef _MSC_VER
#pragma warning(disable:4127)
#endif

// DEBUGON macro variable exists to allow compiler parsing of
// debugprint calls even when DEBUG is false

#if defined(DEBUG) && !defined(NODEBUG)
#define DEBUGON 1
#else
#define DEBUGON 0
#endif

// Define fmt, define file:line:
#define dbgp_prepare_(fmt, xs)                                  \
 char dbgp_fmt_[] = {"%" #fmt};                                 \
 char dbgp_fili_[99];                                           \
 size_t dbgp_nldr_,dbgp_nvec_=5;                                \
 int dbgp_blk_;                                                 \
 snprintf(dbgp_fili_, 99, "%s:%d: ", __FILE__, __LINE__);       \
 dbgp_nldr_ = strlen(xs) + 1;                                   \
 (void) dbgp_nldr_;                                             \
 (void) dbgp_fmt_;                                              \
 (void) dbgp_nvec_;                                             \
 if (sscanf(#fmt,"%d",&dbgp_blk_) <= 0) dbgp_blk_ = 1;          \
 printf("%s", dbgp_fili_);

#define dbgp_conclude_()                        \
 printf("\n");                                  \
 fflush(stdout);
    
// Print only trace info                        
#define print0()                                \
 do { if (DEBUGON) {                            \
     dbgp_prepare_(fmt,"");                     \
     dbgp_conclude_();                          \
   }} while(0)

// Print header                                 
#define printh(s)                               \
 do {                                           \
   if (DEBUGON) {                               \
     dbgp_prepare_(fmt,"");                     \
     printf("%s", s);                           \
     dbgp_conclude_();                          \
   }                                            \
 } while(0)

// Print string                                 
#define prints(x)                               \
 do { if (DEBUGON) {                            \
     dbgp_prepare_(fmt,"");                     \
     printf("%s", #x "=");                      \
     printf("%s", x);                           \
     dbgp_conclude_();                          \
   }} while(0)

// Print variable (e.g. print(d,i) or print(.4f,x))
#define print(fmt,x)                            \
 do {                                           \
   if (DEBUGON) {                               \
     dbgp_prepare_(fmt,#x);                     \
     printf("%s", #x "=");                      \
     printf(dbgp_fmt_, x);                      \
     dbgp_conclude_();                          \
   }                                            \
 } while(0)

// Print 2 variables (e.g. print2(d,i,j))
#define print2(fmt,x,y)                                 \
 do { if (DEBUGON) {                                    \
     dbgp_prepare_(fmt,#x);                             \
     printf("%s", #x "=");                              \
     printf(dbgp_fmt_, (x));                            \
     printf(", " #y "="); printf(dbgp_fmt_, (y));       \
     dbgp_conclude_();                                  \
   }} while(0)

// Print 3 variables (e.g. print3(d,i,j,k))
#define print3(fmt,x,y,z)                       \
 do { if (DEBUGON) {                            \
     dbgp_prepare_(fmt,#x);                     \
     printf("%s", #x "=");                      \
     printf(dbgp_fmt_, x);                      \
     printf(", " #y "="); printf(dbgp_fmt_, y); \
     printf(", " #z "="); printf(dbgp_fmt_, z); \
     dbgp_conclude_();                          \
   }} while(0)

// Print 4 variables
#define print4(fmt,x,y,z,u)                     \
 do { if (DEBUGON) {                            \
     dbgp_prepare_(fmt,#x);                     \
     printf("%s", #x "=");                      \
     printf(dbgp_fmt_, x);                      \
     printf(", " #y "="); printf(dbgp_fmt_, y); \
     printf(", " #z "="); printf(dbgp_fmt_, z); \
     printf(", " #u "="); printf(dbgp_fmt_, u); \
     dbgp_conclude_();                          \
   }} while(0)

// Print 5 variables
#define print5(fmt,x,y,z,u,v)                   \
 do { if (DEBUGON) {                            \
     dbgp_prepare_(fmt,#x);                     \
     printf("%s", #x "=");                      \
     printf(dbgp_fmt_, x);                      \
     printf(", " #y "="); printf(dbgp_fmt_, y); \
     printf(", " #z "="); printf(dbgp_fmt_, z); \
     printf(", " #u "="); printf(dbgp_fmt_, u); \
     printf(", " #v "="); printf(dbgp_fmt_, v); \
     dbgp_conclude_();                          \
   }} while(0)

// Print 6 variables
#define print6(fmt,x,y,z,u,v,w)                 \
 do { if (DEBUGON) {                            \
     dbgp_prepare_(fmt,#x);                     \
     printf("%s", #x "=");                      \
     printf(dbgp_fmt_, x);                      \
     printf(", " #y "="); printf(dbgp_fmt_, y); \
     printf(", " #z "="); printf(dbgp_fmt_, z); \
     printf(", " #u "="); printf(dbgp_fmt_, u); \
     printf(", " #v "="); printf(dbgp_fmt_, v); \
     printf(", " #w "="); printf(dbgp_fmt_, w); \
     dbgp_conclude_();                          \
   }} while(0)

#define printnl()                               \
 if (DEBUGON) printf("\n")

// Print n element vector, dbgp_nvec_ elements per line
// (printv(d,x,7) might print: x=2 3 5 7 11
//                                13 17
#define printv(fmt,x,n)                                                         \
 do { if (DEBUGON) {                                                            \
     dbgp_prepare_(fmt,#x);                                                     \
     printf("%s", #x "=");                                                      \
     int dbgp_idxi_;                                                            \
     for(dbgp_idxi_=0; dbgp_idxi_<(int)(n); dbgp_idxi_++) {                     \
       if (dbgp_idxi_>0 && dbgp_idxi_%dbgp_nvec_==0)                            \
         printf("%s%*s", dbgp_fili_, (int)dbgp_nldr_, "");                      \
       printf(dbgp_fmt_, *((x)+dbgp_idxi_));                                    \
       if (dbgp_idxi_<(int)(n)-1 && dbgp_idxi_%dbgp_nvec_==dbgp_nvec_-1)        \
         printf("\n");                                                          \
       else if (dbgp_idxi_%dbgp_nvec_<dbgp_nvec_-1)                             \
         printf(", ");                                                          \
     }                                                                          \
     dbgp_conclude_();                                                          \
   }} while(0)

// Print m by n matrix in column major order
// printm(5.1f,A,2,2) might print: A= 14.1 -45.2
//                                   111.6  12.8
#define printm(fmt,x,m,n,ix)                                            \
 do { if (DEBUGON) {                                                    \
     dbgp_prepare_(fmt,#x);                                             \
     printf("%s", #x "=");                                              \
     int dbgp_idxi_, dbgp_idxj_;                                        \
     if (n)                                                             \
       for(dbgp_idxi_=0; dbgp_idxi_<(int)(m); dbgp_idxi_++) {           \
         if (dbgp_idxi_>0)                                              \
           printf("%s%*s", dbgp_fili_, (int)dbgp_nldr_, "");            \
         for(dbgp_idxj_=0; dbgp_idxj_<(int)(n); dbgp_idxj_++) {         \
           printf(dbgp_fmt_, *((x)+dbgp_idxj_*(ix)+dbgp_idxi_));        \
           if (dbgp_idxi_<(int)(m)-1 && dbgp_idxj_%(n)==(n)-1)          \
             printf("\n");                                              \
           else                                                         \
             printf(" ");                                               \
         }                                                              \
       }                                                                \
     dbgp_conclude_();                                                  \
   }} while(0)

// print lower triangular matrix
#define printmL(fmt,x,m,n,ix)                                                   \
 do { if (DEBUGON) {                                                            \
     dbgp_prepare_(fmt,#x);                                                     \
     printf("%s", #x "=");                                                      \
     int dbgp_idxi_, dbgp_idxj_;                                                \
     if (n)                                                                     \
       for(dbgp_idxi_=0; dbgp_idxi_<(int)(m); dbgp_idxi_++) {                   \
         if (dbgp_idxi_>0)                                                      \
           printf("%s%*s", dbgp_fili_, (int)dbgp_nldr_, "");                    \
         for(dbgp_idxj_=0; dbgp_idxj_<=(int)(dbgp_idxi_); dbgp_idxj_++) {       \
           printf(dbgp_fmt_, *((x)+dbgp_idxj_*(ix)+dbgp_idxi_));                \
           if (dbgp_idxi_<(int)(m)-1 && dbgp_idxj_ == dbgp_idxi_)               \
             printf("\n");                                                      \
           else                                                                 \
             printf(" ");                                                       \
         }                                                                      \
       }                                                                        \
     dbgp_conclude_();                                                          \
   }} while(0)

// print upper triangular matrix
#define printmU(fmt,x,m,n,ix)                                           \
 do { if (DEBUGON) {                                                    \
     dbgp_prepare_(fmt,#x);                                             \
     printf("%s", #x "=");                                              \
     int dbgp_idxi_, dbgp_idxj_;                                        \
     for(dbgp_idxi_=0; dbgp_idxi_<(int)(m); dbgp_idxi_++) {             \
       if (dbgp_idxi_>0)                                                \
         printf("%s%*s", dbgp_fili_, (int)dbgp_nldr_, "");              \
       for(dbgp_idxj_=0; dbgp_idxj_<dbgp_idxi_; dbgp_idxj_++) {         \
         printf("%*s", dbgp_blk_+1, "");                                \
       }                                                                \
       for(dbgp_idxj_=dbgp_idxi_; dbgp_idxj_<(int)(n); dbgp_idxj_++) {  \
         printf(dbgp_fmt_, *((x)+dbgp_idxj_*(ix)+dbgp_idxi_));          \
         if (dbgp_idxi_<(int)(m)-1 && dbgp_idxj_%(n)==(n)-1)            \
           printf("\n");                                                \
         else                                                           \
           printf(" ");                                                 \
       }                                                                \
     }                                                                  \
     dbgp_conclude_();                                                  \
   }} while(0)

// Print m by n matrix in row major order
#define printmR(fmt,x,m,n,ix)                                   \
 do { if (DEBUGON) {                                            \
     dbgp_prepare_(fmt,#x);                                     \
     printf("%s", #x "=");                                      \
     int dbgp_idxi_, dbgp_idxj_;                                \
     for(dbgp_idxi_=0; dbgp_idxi_<(int)(m); dbgp_idxi_++) {     \
       if (dbgp_idxi_>0)                                        \
         printf("%s%*s", dbgp_fili_, (int)dbgp_nldr_, "");      \
       for(dbgp_idxj_=0; dbgp_idxj_<(int)(n); dbgp_idxj_++) {   \
         printf(dbgp_fmt_, *((x)+dbgp_idxi_*(ix)+dbgp_idxj_));  \
         if (dbgp_idxi_<(int)(m)-1 && dbgp_idxj_%(n)==(n)-1)    \
           printf("\n");                                        \
         else                                                   \
           printf(" ");                                         \
       }                                                        \
     }                                                          \
     dbgp_conclude_();                                          \
   }} while(0)

// Simpler print functions (here fmt must be enclosed in "")
#define debugprint(fmt, ...)                                    \
 do { if (DEBUGON) {                                            \
     printf("%s:%d: " fmt, __FILE__, __LINE__, __VA_ARGS__);    \
   }} while (0)

#define dp_leader()                             \
 do { if (DEBUGON) {                            \
     printf("%s:%d: ", __FILE__, __LINE__);     \
   }} while (0)

#define dp_print(fmt, ...)                      \
 do { if (DEBUGON) {                            \
     printf(fmt, __VA_ARGS__);                  \
   }} while (0)

#define dp_string(s)                            \
 do { if (DEBUGON) {                            \
     printf("%s",s);                            \
   }} while (0)

#define dp_newline()                            \
 do { if (DEBUGON) {                            \
     printf("\n");                              \
   }} while (0)

#endif

