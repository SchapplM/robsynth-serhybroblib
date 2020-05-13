% Geometrische Jacobi-Matrix für beliebiges Segment von
% hybBKplanar
% Use Code from Maple symbolic Code Generation
% 
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,CF,ED]';
% 
% Output:
% Jg [6x3]
%   Geometrische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 19:03
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = hybBKplanar_jacobig_sym_varpar(qJ, link_index, r_i_i_C, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'hybBKplanar_jacobig_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'hybBKplanar_jacobig_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'hybBKplanar_jacobig_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'hybBKplanar_jacobig_sym_varpar: pkin has to be [6x1] (double)');

% Function calls
Ja_transl = hybBKplanar_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin);
Jg_rot = hybBKplanar_jacobig_rot_sym_varpar(qJ, link_index, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
end