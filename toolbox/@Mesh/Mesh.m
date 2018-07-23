classdef Mesh < handle
%MESH Summary of this function goes here
%   Detailed explanation goes here

properties
    F
    nF
    Nf
    V
    nV
    Nv
    F2V
    V2E
    A
    E
    nE
    BE
    BV
    Aux
    V2V
    E2E
    E2F
end

methods
    function obj = Mesh(varargin)
        if (length(varargin)==1) && isa(varargin{1},'Mesh')
            obj.V=varargin{1}.V;
            obj.nV=varargin{1}.nV;
            obj.Nv=varargin{1}.Nv;
            obj.F=varargin{1}.F;
            obj.nF=varargin{1}.nF;
            obj.Nf=varargin{1}.Nf;
            obj.F2V=varargin{1}.F2V;
            obj.V2E=varargin{1}.V2E;
            obj.A=varargin{1}.A;
            obj.E=varargin{1}.E;
            obj.nE=varargin{1}.nE;
            obj.BE=varargin{1}.BE;
            obj.BV=varargin{1}.BV;
            obj.Aux=varargin{1}.Aux;
            obj.V2V=varargin{1}.V2V;
            obj.E2E=varargin{1}.E2E;
            obj.E2F=varargin{1}.E2F;
        elseif length(varargin)>=2
            obj.Aux.filename = varargin{2};
            switch(varargin{1})
                case 'off'
                    [obj.V,obj.F] = obj.read_off(varargin{2});
                case 'obj'
                    [obj.V,obj.F] = obj.read_obj(varargin{2});
                case 'VF'
                    if size(varargin{2},1)==2
                        obj.V=[varargin{2};zeros(1,size(varargin{2},2))];
                    else
                        obj.V=varargin{2};
                    end
                    obj.F=varargin{3};
            end
            obj.nV = size(obj.V,2);
            obj.nF = size(obj.F,2);
            obj.F2V = obj.ComputeF2V();
            obj.V2E = obj.ComputeV2E();
            obj.nE = size(obj.V2E,2);
            obj.E = obj.ComputeE();
            obj.E2F = obj.ComputeE2F();
            
        else
            obj.F=[];
            obj.V=[];
        end
    end
end
methods(Static)
    [V,F,Fs] = read_off(filename)
    [V,F,Fs] = read_obj(filename)
    v = getoptions(options, name, v, mendatory)
    [D,S,Q] = PerformFrontPropagation(vertex, faces, W,start_points,end_points, nb_iter_max, H, L, values, dmax)
end
end


