classdef Rng < handle
  %VARMAPACK.RNG  Randompack generator for Varmapack simulation.
  %
  %   rng = VARMAPACK.RNG() creates a randomized generator using the
  %   default Randompack engine. rng = VARMAPACK.RNG(seed) creates the
  %   same deterministic stream whenever seed is reused.

  properties (Access = private)
    handle = uint64(0)
  end

  methods
    function obj = Rng(seed)
      obj.handle = varmapack_rng_gateway('create');
      if nargin > 0
        try
          obj.seed(seed);
        catch exception
          delete(obj);
          rethrow(exception);
        end
      end
    end

    function seed(obj, seed)
      %SEED  Reset this generator to a deterministic stream.
      varmapack_rng_gateway('seed', obj.mex_handle(), seed);
    end

    function delete(obj)
      if obj.handle ~= 0
        varmapack_rng_gateway('free', obj.handle);
        obj.handle = uint64(0);
      end
    end
  end

  methods (Hidden)
    function handle = mex_handle(obj)
      if ~isscalar(obj) || obj.handle == 0
        error('varmapack:Rng:invalid', 'RNG object is no longer valid');
      end
      handle = obj.handle;
    end
  end
end
