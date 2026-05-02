macOS: Uses arc4random_buf(), which Apple documents as cryptographically secure. This approach is used by Go's crypto/rand package and is a fallback in Python's os.urandom(). While some libraries like Rust's getrandom (using getentropy()) or libsodium (using /dev/urandom for maximum conservatism), arc4random_buf() is a native, well-supported system API suitable for cryptographic purposes.

Windows: Uses RtlGenRandom (exported from advapi32.dll as SystemFunction036). This function is used by Microsoft's own C runtime rand_s() and by major projects including libsodium, Rust, Firefox, Chrome/Chromium, and BoringSSL. On modern Windows (Vista+), it internally calls ProcessPrng which uses AES-CTR-DRBG as specified by NIST SP800-90. While officially undocumented, it cannot be removed without breaking Microsoft's C standard library implementation. Uses the algorithm AES-CTR 

SHORT: Windows: Uses RtlGenRandom (advapi32.dll, SystemFunction036). Used by libsodium, Rust, Firefox, Chrome, and Microsoft's own rand_s(). Internally calls ProcessPrng on modern Windows (AES-CTR-DRBG, NIST SP800-90 compliant).

OpenBSD: Uses arc4random_buf(). OpenBSD's implementation uses ChaCha20 (since OpenBSD 5.5, 2014) and is documented by OpenBSD as the recommended interface for cryptographic randomness. Used by Go crypto/rand and explicitly trusted by libsodium (marked as HAVE_SAFE_ARC4RANDOM). OpenBSD's documentation recommends arc4random(3) over getentropy(2) for general use.

FreeBSD: Uses getrandom() syscall (available since FreeBSD 12.0, 2018). Used by Go crypto/rand. Modern FreeBSD uses ChaCha20 since FreeBSD 12.0.

Linux: Uses the getrandom() system call (available since Linux 3.17, 2014) when available, otherwise falls back to /dev/urandom. Both use the kernel's ChaCha20-based CSPRNG (since kernel 5.6, 2020; earlier kernels used other algorithms). The implementation follows the same pattern as libsodium, using a syscall to avoid an undefined symbol link error, viz:

#if defined(__linux__) && defined(SYS_getrandom)
  ssize_t result = syscall(SYS_getrandom, buffer, size, 0);
  if (result != size) {
    // Kernel doesn't support it, fall back to /dev/urandom
  }

SHORT: Linux: Uses getrandom() syscall when available (Linux 3.17+), with /dev/urandom fallback for older kernels. Uses ChaCha20-based CSPRNG (kernel 5.6+).

All these systems implement reseeding.

Here is the platform detection skeleton:
#if defined(_WIN32)
  // Windows: RtlGenRandom
#elif defined(__APPLE__)
  // macOS/iOS: arc4random_buf
#elif defined(__OpenBSD__)
  // OpenBSD: arc4random_buf
#elif (defined(__linux__) || defined(__FreeBSD__) || defined(__DragonFly__)) && defined(SYS_getrandom)
  // Linux/FreeBSD/DragonFly: getrandom syscall
#else
  // Fallback: /dev/urandom (works on Solaris, illumos, NetBSD, old systems)
#endif

Our implementation is primarily drawn from libsodium with additional reference to Go (crypto/rand), Rust (getrandom crate), Python (os.urandom), Firefox/Mozilla, Chrome/Chromium.



xoshiro’s jump functions aren’t exposed because splitmix64 seeding already gives statistically independent 256-bit starting states with essentially zero collision probability. Melissa O’Neill has noted the same principle for PCG (independent seeding rather than jump APIs), and NumPy likewise chose not to expose xoshiro’s jump and long-jump operations. 


PCG64’s 128-bit state gives a cycle of 3.4×10³⁸, so even extreme HPC workloads consume an utterly negligible fraction of it — effectively guaranteeing no overlaps or repeats in any real simulation.
Different increments theoretically give mathematically disjoint orbits, but in practice distinct seeds alone already push streams so far apart that using a fixed increment is completely safe, as also done in NumPy’s PCG64.
Lemire style = the trick to avoid bias in random integers.
Lemire's algorithm:
cuint32_t random_bounded(uint32_t n) {
  uint64_t m = (uint64_t)random_uint32() * (uint64_t)n;
  uint32_t l = (uint32_t)m;
  
  if (l < n) {
    uint32_t t = -n % n;  // or: (UINT32_MAX - n + 1) % n
    while (l < t) {
      m = (uint64_t)random_uint32() * (uint64_t)n;
      l = (uint32_t)m;
    }
  }
  
  return m >> 32;  // High 32 bits
}
How it works:

Multiply random 32-bit value by n → gives 64-bit result
The high 32 bits give the answer (usually)
Rejection sampling fixes the bias when needed (rare)

Benefits:

Fast - multiplication instead of division
Unbiased - mathematically correct
Simple - no lookup tables

Published by: Daniel Lemire in 2019: "Fast Random Integer Generation in an Interval"
This is what modern libraries like PCG and many game engines use. It's much faster than randombytes_uniform() style rejection sampling with modulo.
What I'm thinking of adding is: 
- fairly comprehensive tests
- timing functions
- float support
- more distributions (exponential, gamma, beta,...)
- features to extract – restore the rng (we discussed that at some point)

One thing I want to investigate now: I'd like to time (and even compare) ziggurat normals with polar normals. I already have the polar normal function. How would I design an interface to offer it?
