function DrawDemo
  draw_u01('x256++', 42, 2);
  draw_u01('pcg64', 42, 2);
  draw_philox(1, 42);
  draw_u01('chacha20', 42, 2);
  draw_squares(1, 123456789);
end

function draw_u01(engine, seed, key)
  rng = randompack_create(engine);
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_seed(rng, seed, uint32(key));
  x = randompack_u01(rng, 2);
  fprintf("%s seed=%d spawnkey=%d u01=%.17g %.17g\n", engine, seed, key, x(1), x(2));
end

function draw_philox(ctr0, key0)
  rng = randompack_create('philox');
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_set_state(rng, uint64([ctr0 0 0 0 0 0]));
  randompack_philox_set_key(rng, uint64([key0 0]));
  x = randompack_u01(rng, 2);
  fprintf("%s counter=%d key=%d u01=%.17g %.17g\n", 'philox', ctr0, key0, x(1), x(2));
end

function draw_squares(ctr, key)
  rng = randompack_create('squares');
  cleanup = onCleanup(@() randompack_free(rng));
  randompack_set_state(rng, uint64([ctr 0]));
  randompack_squares_set_key(rng, key);
  x = randompack_u01(rng, 2);
  fprintf("%s counter=%d key=%d u01=%.17g %.17g\n", 'squares', ctr, key, x(1), x(2));
end
