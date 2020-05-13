% Analytische Jacobi-Matrix für beliebiges Segment von
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in picker2Dm1OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% Ja [6x12]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = picker2Dm1OL_jacobia_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_jacobia_sym_varpar: qJ has to be [12x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'picker2Dm1OL_jacobia_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'picker2Dm1OL_jacobia_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_jacobia_sym_varpar: pkin has to be [8x1] (double)');

% Function calls
Ja_transl = picker2Dm1OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin);
Ja_rot = picker2Dm1OL_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
end