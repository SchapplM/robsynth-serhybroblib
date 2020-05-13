% Geometrische Jacobi-Matrix für beliebiges Segment von
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% Jg [6x2]
%   Geometrische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = fourbar1turnDE2_jacobig_sym_varpar(qJ, link_index, r_i_i_C, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_jacobig_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'fourbar1turnDE2_jacobig_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE2_jacobig_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_jacobig_sym_varpar: pkin has to be [5x1] (double)');

% Function calls
Ja_transl = fourbar1turnDE2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin);
Jg_rot = fourbar1turnDE2_jacobig_rot_sym_varpar(qJ, link_index, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
end