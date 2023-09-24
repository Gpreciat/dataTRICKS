% Support only 2 dim array
classdef qdouble
   properties (Constant)
      mex = qdouble.mexSelector();
   end
   
   methods (Static)
      function func = mexSelector()
         if isarm()
            func = @qdoubleArmMex;
         else
            func = @qdoubleMex;
         end
      end
   end
   
   properties
      x;
   end
   
   methods
      function o = qdouble(a)
         if iscell(a) || isa(a, 'uint8')
            o.x = a;
         else
            o.x = qdouble.toMex(a);
         end
      end
      
      %% Sizes
      function [m, n] = size(a, dim)
         a = a.x;
         if (iscell(a))
            s = double(size(a{1}));
         else
            s = [size(a, 2) size(a, 3)];
         end
         
         if (nargin == 1)
            if (nargout == 2)
               m = s(1);
               n = s(2);
            else
               m = s;
            end
         else
            m = s(dim);
         end
      end
      
      function r = length(a)
         s = size(a);
         if (min(s) == 0)
            r = 0;
         else
            r = max(s);
         end
      end
      
      function r = numel(a)
         r = prod(size(a)); %#ok<*PSIZE>
      end
      
      function r = isscalar(a)
         r = (numel(a) == 1);
      end
      
      function r = isvector(a)
         r = (size(a,1) == 1 || size(a,2) == 1);
      end
      
      function r = isempty(a)
         r = (numel(a) == 0);
      end
      
      function r = isrow(a)
         r = (size(a,1) == 1);
      end
      
      function r = iscolumn(a)
         r = (size(a,2) == 1);
      end
      
      function r = ismatrix(~)
         r = true;
      end
      
      function r = end(a, k, n)
         if n == 2
            r = size(a,k);
         else
            r = size(a,1) * size(a,2);
         end
      end
      
      %% Conversions
      function disp(a)
         disp(double(a))
      end
      
      function r = double(a)
         r = qdouble.UnaryOp('double', a);
      end
      
      function r = logical(a)
         r = qdouble.UnaryOp('logical', a);
      end
      
      function r = sparse(a,b,c,d,e,f)
         if nargin == 1
            if ~issparse(a)
               A = logical(sparse(double(a)));
               s = struct('type', '()');
               s.subs = {find(A)};
               s = qdouble(subsref(a, s));
               s.x = reshape(s.x, [size(s.x,1) numel(s.x)/size(s.x,1)]);
               r = qdouble({A;s.x});
            else
               r = qdouble(a);
            end
         else
            assert(nargin == 3 || nargin == 5 || nargin == 6, "unsupported number of arguments for sparse");
            if nargin == 3
               i = a; j = b; s = c; m = max(i); n = max(j); nzmax = length(s);
               if isempty(i)
                  m = 0; n = 0;
               end
            elseif nargin == 5
               i = a; j = b; s = c; m = d; n = e; nzmax = length(s);
            elseif nargin == 6
               i = a; j = b; s = c; m = d; n = e; nzmax = f;
            end
            A = sparse(i, j, ones(nzmax,1,'logical'),m,n,nzmax);
            assert(length(s) == nzmax);
            s = qdouble(s(:));
            r = qdouble({A;s.x});
         end
      end
      
      function r = full(a)
         if issparse(a)
            A = a.x{1};
            s = a.x{2};
            r = zeros(size(s,1), size(A,1), size(A,2), 'uint8');
            idx = find(a.x{1});
            r(:, idx) = s;
            r = qdouble(r);
         else
            r = qdouble(a);
         end
      end
      
      function r = issparse(a)
         r = iscell(a.x);
      end
      
      function r = issymmetric(a)
         if size(a,1) ~= size(a,2)
            r = false;
         else
            r = full(all(a == a', 'all'));
         end
      end
      
      function r = diag(a)
         assert(nargin == 1, "Support only one parameters")
         if isvector(a)
            r = qdouble(sparse(a));
            if size(a,1) == 1
               r = r';
            end
            r.x{1} = diag(double(a));
            if ~issparse(a)
               r = full(r);
            end
         else
            k = min(size(a));
            if k == 0
               r = qdouble([]);
            else
               idx = (1:k) + size(a,1) * (0:(k-1));
               r =  subsref(a, struct('type', '()', 'subs', {{idx}}));
            end
            if issparse(a)
               r = sparse(r);
            end
            if size(a,2) < size(a,1)
               r = r';
            end
         end
      end
      
      function r = reshape(a, s)
         if issparse(a)
            a.x{1} = reshape(a.x{1}, s);
         else
            a.x = reshape(a.x, [size(a.x,1) s]);
         end
         r = qdouble(a.x);
      end
      
      function r = eps(~)
         r = qdouble(qdouble.mex('eps'));
      end
      
      %% Comparisons
      function r = lt(a, b)
         r = qdouble.BinaryOp('lt', a, b);
         if ((issparse(a) || issparse(b)) && ~issparse(r)), r = sparse(r); end
      end
      
      function r = gt(a, b)
         r = qdouble.BinaryOp('gt', a, b);
         if ((issparse(a) || issparse(b)) && ~issparse(r)), r = sparse(r); end
      end
      
      function r = le(a, b)
         r = ~gt(a,b);
      end
      
      function r = ge(a, b)
         r = ~lt(a,b);
      end
      
      function r = ne(a, b)
         r = qdouble.BinaryOp('ne', a, b);
         if ((issparse(a) || issparse(b)) && ~issparse(r)), r = sparse(r); end
      end
      
      function r = eq(a, b)
         r = ~ne(a,b);
      end
      
      function r = and(a, b)
         r = qdouble.BinaryOp('and', a, b);
      end
      
      function r = or(a, b)
         r = qdouble.BinaryOp('or', a, b);
      end
      
      function r = not(a)
         r = ~logical(a);
      end
      
      %% Arithmetic
      function r = plus(a, b)
         r = qdouble.BinaryOp('plus', a, b);
      end
      
      function r = minus(a, b)
         r = qdouble.BinaryOp('minus', a, b);
      end
      
      function r = times(a, b)
         r = qdouble.BinaryOp('times', a, b);
      end
      
      function r = power(a, b)
         assert(isscalar(b) && b == 2, "support only power(a, 2)");
         r = a .* a;
      end
      
      function r = rdivide(a, b)
         if issparse(b)
            b = full(b);
         end
         r = qdouble.BinaryOp('rdivide', a, b);
         
         if issparse(a) % To fix the bug regarding x / 0
            r2 = double(a) ./ double(b);
            r(~isfinite(r2)) = r2(~isfinite(r2));
         end
      end
      
      function r = ldivide(a, b)
         r = b./a;
      end
      
      function r = uminus(a)
         r = qdouble.UnaryOp('uminus', a);
      end
      
      function r = uplus(a)
         r = a;
      end
      
      function r = abs(a)
         r = qdouble.UnaryOp('abs', a);
      end
      
      function r = sqrt(a)
         r = qdouble.UnaryOp('sqrt', a);
      end
      
      function r = sum(varargin)
         r = qdouble.ReductionOp('sum', varargin{:});
      end
      
      function r = prod(varargin)
         r = qdouble.ReductionOp('prod', varargin{:});
      end
      
      function r = all(a, dim)
         if nargin == 1
            r = all(logical(a));
         else
            r = all(logical(a), dim);
         end
      end
      
      function r = any(a, dim)
         if nargin == 1
            r = any(logical(a));
         else
            r = any(logical(a), dim);
         end
      end
      
      function r = max(a, b, c)
         if nargin == 1
            r = qdouble.ReductionOp('max', a);
         elseif nargin == 3 && size(b,1) == 0 && size(b,2) == 0
            r = qdouble.ReductionOp('max', a, c);
         else
            assert(nargin == 2);
            r = qdouble.BinaryOp('max2', a, b);
            if issparse(a) || issparse(b)
               r = sparse(r);
            end
         end
      end
      
      function r = min(a, b, c)
         if nargin == 1
            r = qdouble.ReductionOp('min', a);
         elseif nargin == 3 && size(b,1) == 0 && size(b,2) == 0
            r = qdouble.ReductionOp('min', a, c);
         else
            assert(nargin == 2);
            r = qdouble.BinaryOp('min2', a, b);
            if issparse(a) || issparse(b)
               r = sparse(r);
            end
         end
      end
      
      function r = norm(a)
         assert(isvector(a), "norm supports only vectors.");
         r = sqrt(full(sum(a.^2)));
      end
      
      function r = nnz(a)
         r = full(sum(logical(a),'all'));
      end
      
      %% Access
      function [i,j,v] = find(a)
         if nargout == 1
            i = find(logical(a));
         elseif nargout == 2
            [i,j] = find(logical(a));
         else
            a = sparse(qdouble(a));
            [i,j] = find(a.x{1});
            v = qdouble(a.x{2});
            if size(a,1) == 1
               v = v';
            end
            if numel(v) == 0
               v = qdouble([]);
            end
         end
      end
      
      %% Matrix operations
      function r = ctranspose(a)
         r = transpose(a);
      end
      
      function r = transpose(a)
         r = qdouble.UnaryOp('transpose', a);
      end
      
      function r = mtimes(a, b)
         if isscalar(a) || isscalar(b)
            r = a .* b;
         else
            r = qdouble.BinaryOp('mtimes', a, b);
         end
      end
      
      function r = mldivide(a, b)
         assert(size(a,1) == size(b,1), 'Incompatible Size.');
         s = [size(a,2), size(b,2)];
         if size(a,2) == 0 || size(a,1) == 0 || size(b,2) == 0
            r = qdouble(zeros(s));
         elseif isscalar(a)
            r = reshape(b ./ a, s);
         else
            b_full = full(b); % sparsity usually does not help. Do this for simplicity
            r = qdouble.BinaryOp('mldivide', a, b_full);
         end
         if (issparse(a) && issparse(b))
            r = sparse(r);
         end
      end
      
      function r = mrdivide(a, b)
         r = (b'\a')';
      end
      
      function r = horzcat(varargin)
         allDense = true;
         for k = 1:length(varargin)
            allDense = allDense & ~issparse(varargin{k});
         end
         
         if (allDense)
            cat_input = cell(length(varargin), 1);
            for k = 1:length(varargin)
               cat_input{k} = varargin{k}.x;
            end
            r = qdouble(cat(3, cat_input{:}));
         else
            A_bool = [];
            m = []; v = [];
            for k = 1:length(varargin)
               A = sparse(qdouble(varargin{k}));
               if k == 1
                  m = size(A,1);
                  A_bool = A.x{1};
                  v = A.x{2};
               else
                  assert(m == size(A,1), 'Incompatible size.');
                  A_bool = [A_bool A.x{1}];
                  v = [v A.x{2}];
               end
            end
            r = qdouble({A_bool;v});
         end
      end
      
      function r = vertcat(varargin)
         for i = 1:length(varargin)
            varargin{i} = varargin{i}';
         end
         r = horzcat(varargin{:})';
      end
      
      function r = subsref(a,s)
         if strcmp(s(1).type, '()')
            if issparse(a)
               idx = subsref(qdouble.ToIndex(a),s);
               A = logical(idx);
               [~,~,idx] = find(idx);
               v = a.x{2};
               v = v(:, idx);
               r = qdouble({A;v});
            else
               idx = subsref(qdouble.ToIndex(a),s);
               r = qdouble(reshape(a.x(:,idx), [size(a.x,1) size(idx)]));
            end
         else
            r = builtin('subsref',a,s);
         end
      end
      
      function r = subsasgn(a,s,b)
         if strcmp(s(1).type, '()')
            [a_idx, a_v] = qdouble.ToIndex(a);
            [b_idx, b_v] = qdouble.ToIndex(b);
            b_idx = max(a_idx,[],'all') + b_idx;
            v = [a_v b_v];
            
            r_idx = subsasgn(a_idx,s,b_idx);
            
            if issparse(r_idx)
               A = logical(r_idx);
               [~,~,idx] = find(r_idx);
               v = v(:, idx);
               r = {A;v};
            else
               r = reshape(v(:, r_idx), [size(v,1) size(r_idx)]);
            end
            r = qdouble(r);
         else
            r = builtin('subsasgn',a,s,b);
         end
      end
      
      %% cholesky related
      function r = dissect(a)
         r = dissect(double(a));
      end
      
      function r = amd(a)
         r = amd(double(a));
      end
      
      function r = colamd(a)
         r = colamd(double(a));
      end
      
      function r = symbfact(a)
         r = symbfact(double(a));
      end
      
      function [r1, r2] = etree(a)
         [r1, r2] = etree(double(a));
      end
      
      function r = chol(a)
         r = qdouble.UnaryOp('chol', a);
      end
   end
   
   methods (Static)
      function r = ones(varargin)
         r = qdouble(ones(varargin{:}));
      end
      
      function r = zeros(varargin)
         r = qdouble(zeros(varargin{:}));
      end
      
      function r = eye(varargin)
         r = qdouble(eye(varargin{:}));
      end
      
      function r = rand(varargin)
         r = qdouble(rand(varargin{:}));
      end
      
      function r = randn(varargin)
         r = qdouble(randn(varargin{:}));
      end
      
      function r = randi(varargin)
         r = qdouble(randi(varargin{:}));
      end
      
      function r = toMex(a)
         if isa(a, 'qdouble')
            r = a.x;
         elseif isobject(a)
            r = a.x;
            r = qdouble.mex('toMex', r);
         else
            r = qdouble.mex('toMex', double(a));
         end
      end
      
      function r = UnaryOp(cmd, a)
         a = qdouble.toMex(a);
         r = qdouble.mex(cmd, a);
         
         if iscell(r) || isa(r, 'uint8')
            r = qdouble(r);
         end
      end
      
      function r = BinaryOp(cmd, a, b)
         a = qdouble.toMex(a);
         b = qdouble.toMex(b);
         r = qdouble.mex(cmd, a, b);
         
         if iscell(r) || isa(r, 'uint8')
            r = qdouble(r);
         end
      end
      
      function [r, v] = ToIndex(a)
         if nargout == 2 && ~isa(a, 'qdouble')
            a = qdouble(a);
         end
         
         if issparse(a)
            [i,j] = find(double(a));
            r = sparse(i,j,1:length(i),size(a,1),size(a,2));
            
            if nargout == 2
               v = a.x{2};
            end
         else
            r = reshape(1:numel(a),[size(a,1) size(a,2)]);
            
            if nargout == 2
               v = reshape(a.x, [size(a.x,1) numel(a.x)/size(a.x,1)]);
            end
         end
      end
      
      function r = ReductionOp(cmd, a, dim_)
         a_size = size(a);
         a = qdouble(a);
         
         if nargin == 2
            if (size(a,1) > 1)
               dim = 1;
            elseif (size(a,2) > 1)
               dim = 2;
            else
               dim = 1;
               a = subsref(a, struct('type', '()', 'subs', {{':'}}));
            end
         else
            assert(strcmp(dim_, 'all') || isscalar(dim_), 'dim can either be all or a scalar');
            if strcmp(dim_, 'all')
               dim = 1;
               a = subsref(a, struct('type', '()', 'subs', {{':'}}));
            else
               dim = dim_;
            end
         end
         
         if size(a,1) == 0
            switch cmd
               case 'sum'
                  r = qdouble(0.0 * ones(1, size(a,2)));
               case 'prod'
                  r = qdouble(1.0 * ones(1, size(a,2)));
               case 'max'
                  r = qdouble(ones(0, a_size(2)));
               case 'min'
                  r = qdouble(ones(0, a_size(2)));
            end
            if dim == 2
               r = r';
            end
         elseif dim == 1
            r = qdouble.mex(cmd, a.x);
            r = qdouble(r);
         elseif dim == 2
            a = a';
            r = qdouble.mex(cmd, a.x);
            r = qdouble(r)';
         end
         
         if issparse(a)
            r = sparse(r);
         end
      end
   end
end