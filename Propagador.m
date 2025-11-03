classdef Propagador < handle
    properties
        prop
        N
        nz

        counter = 0;
    end

    methods
        function obj = Propagador(prop, N, nz, counter)
            if nargin > 0
                obj.prop = prop;
                obj.N = N;
                obj.nz = nz;
            end
        end
        
        function Uz = propagar(obj, U0)
            obj.counter = obj.counter + 1;
            Uz = zeros(obj.N, obj.N, obj.nz);
            Uz(:,:,1) = U0;
            F = fftshift(fft2(fftshift(U0)));
            for ii=2:obj.nz
                F = F .* obj.prop;
                A = ifftshift(ifft2(ifftshift(F)));
                Uz(:,:,ii) = A;
            end
        end
    end
end