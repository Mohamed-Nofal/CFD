function e=tri_sol(a,b,c,d,M)
% '____________________________________________________________
% tri:     ' Subroutine determine the
% '          solution of a Tridiagonal matrix given by
% '          b(i)*e(i-1) + d(i)*e(i) + a(i)*e(i+1) = c(i)
% '          **********************************************
           
            for i=2:M
            t = b(i) / d(i - 1);
            d(i) = d(i) - t * a(i - 1);
            c(i) = c(i) - t * c(i - 1);
            end
                e(M) = c(M) / d(M);
            for k=2:M
            i = M - k + 1;
            e(i) = (c(i) - a(i) * e(i + 1)) / d(i);
            end
           % RETURN
% '____________________________________________________________________
