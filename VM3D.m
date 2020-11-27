classdef VM3D < handle
  %------------------------------------------------------------------------
  properties
    young;
    poisson;
    SY = TABLE.empty;
    TOLL = 0.0000001;
    MAXITER = 1000;
  end
  %------------------------------------------------------------------------
  methods
    %----------------------------------------------------------------------
    function this = VM3D(young, poisson, yield_stress)
      this.young = young;
      this.poisson = poisson;
      this.SY = yield_stress;
    end
    %----------------------------------------------------------------------
    function [I1] = firstinv(this, var)
      I1 = (var(1) + var(2) + var(3));
    end
    %----------------------------------------------------------------------
    function [J2] = secondinvdev(this, var)
      J2 = (1 / 3) * ((var(1) - var(2))^2 + ...
                      (var(2) - var(3))^2 + ...
                      (var(3) - var(1))^2 + ...
                      6 * (var(4)^2 + var(5)^2 + var(6)^2));
    end
    %----------------------------------------------------------------------
    function [n] = norm(this, vec)
      n = vec' * vec;
      n = sqrt(n);
    end
    %----------------------------------------------------------------------
    function [dJ2] = dsecondinvdev(this, var)
      dJ2 = [+(var(1) - var(2)) - (var(3) - var(1));
             -(var(1) - var(2)) + (var(2) - var(3));
             -(var(2) - var(3)) + (var(3) - var(1));
             6.0 * var(4);
             6.0 * var(5);
             6.0 * var(6);] * (2 / 3);
    end
    %----------------------------------------------------------------------
    function [k] = eqstrain(this, pstrain)
      k = sqrt(this.secondinvdev(pstrain));
    end
    %----------------------------------------------------------------------
    function [dk] = deqstrain(this, pstrain, k, df_ds)
      dk = 0.0;
      if (k == 0.0)
        return
      end
      dJ2 = (this.dsecondinvdev(pstrain));
      dstrain_dg = df_ds;
      dk = (dJ2' * dstrain_dg) / k;
    end
    %----------------------------------------------------------------------
    function [f] = get_f(this, stress, k)
      I1 = this.firstinv(stress);
      N  = this.norm(stress);
      sy = this.SY.getTableVal(k);
      f = sqrt(N^2 - (I1^2) / 3) - sy;
    end
    %----------------------------------------------------------------------
    function [df] = get_dfds(this, stress)
      df = zeros(6, 1);
      I1 = this.firstinv(stress);
      N = this.norm(stress);
      den = sqrt(N^2 - (I1^2) / 3);
      if (den ~= 0)
        I = [1;1;1;0;0;0;];
        df = stress - (I1 / 3) * I;
        df = df / den;
      end
    end
    %----------------------------------------------------------------------
    function [ddf] = get_ddfdds(this, stress)
      ddf = zeros(6, 6);
      I1 = this.firstinv(stress);
      N = this.norm(stress);
      den = N^2 - (I1^2) / 3;
      if (den <= 0)
        return
      end
      II = eye(6, 6);
      I = [1;1;1;0;0;0;];
      ddf = (II - I * I') / sqrt(den);
      ddf = ddf-((stress - I1 * I) * (stress - (I1 / 3) * I)') * den^-1.5;
    end
    %----------------------------------------------------------------------
    function [df] = get_dfdgamma(this, k, dk_dgamma)
      df = 0.0;
      if (k == 0.0)
        return
      end
      df_dsy = -1;
      dsy_dk = this.SY.getTableDVal(k);
      df = df_dsy * dsy_dk * dk_dgamma;
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    function [D] = GETELASTICCONSTMAT(this)
      E  = this.young;
      nu = this.poisson;
      D  = zeros(6, 6);
      v01 = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
      v02 = v01 * (1.d0 - 2.d0 * nu) / 2.0;
      D(1,1) = (1.0-nu)*v01; D(1,2) = v01*nu;       D(1,3) = v01*nu;
      D(2,1) = v01*nu;       D(2,2) = (1.0-nu)*v01; D(2,3) = v01*nu;
      D(3,1) = v01*nu;       D(3,2) = v01*nu;       D(3,3) = (1.0-nu)*v01;
      D(4,4) = v02; 
      D(5,5) = v02; 
      D(6,6) = v02;
    end
    %----------------------------------------------------------------------
    function [stress, Dep, statev] = GETCONSTMAT(this, ...
                                            strain, dstrain, statev)
      check = true;
      strain = strain + dstrain;
      %elastic trial
      pstrain = statev(1:6);
      ko = statev(7);
      k  = ko;
      Dep = this.GETELASTICCONSTMAT();
      stress = Dep * (strain - pstrain);
      f = this.get_f(stress, k);
      %plastic check
      if (f > 0.0)
        II = eye(6, 6);
        pstraino = pstrain;
        trStress = stress;
        gamma = 0.0;
        r = zeros(6, 1);
        res = 1.0;
        x = [stress; gamma;];
        ITER = 0;
        while ((abs(f) > this.TOLL) || (res > this.TOLL))
          h = [f; r;];
          df_ds = this.get_dfds(stress);
          ddf_dds = this.get_ddfdds(stress);
          dk_dg = this.deqstrain(pstrain, k, df_ds);
          df_dg = this.get_dfdgamma(k, dk_dg);
          dr_ds = -II - gamma * Dep * ddf_dds;
          dr_dg = -Dep * df_ds;
          B = [df_ds', df_dg;
               dr_ds, dr_dg;];
          dx = -B \ h;
          x = x + dx;
          stress(1:6) = x(1:6);
          gamma = x(7);
          dpstrain = gamma * df_ds;
          pstrain = pstraino + gamma * df_ds;
          k = ko + this.eqstrain(dpstrain);
          f = this.get_f(stress, k);
          r = trStress - stress - Dep * dpstrain;
          res = this.norm(r);
          ITER = ITER + 1;
          fprintf('  ITER: %i \t\t f: %e \t\t res: %e \n', ITER, f, res);
          if (ITER > this.MAXITER)
            check = false;
            stress = Dep * (strain - dstrain - pstraino);
            fprintf('WARNING: convergence not found\n')
            break
          end
        end
      end
      statev(1:6) = pstrain(1:6);
      statev(7) = k;
      statev(8) = sqrt(this.secondinvdev(stress));
      %total equivalent strain evaluate using the total strain
      statev(9) = this.eqstrain(strain);
      statev(10) = check;
    end
    %----------------------------------------------------------------------
  end
  %------------------------------------------------------------------------
end