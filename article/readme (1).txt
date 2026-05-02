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
