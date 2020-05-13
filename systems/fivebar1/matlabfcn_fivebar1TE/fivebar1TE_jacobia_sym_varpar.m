% Analytische Jacobi-Matrix für beliebiges Segment von
% fivebar1TE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in fivebar1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% Ja [6x2]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = fivebar1TE_jacobia_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_jacobia_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'fivebar1TE_jacobia_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fivebar1TE_jacobia_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_jacobia_sym_varpar: pkin has to be [5x1] (double)');

% Function calls
Ja_transl = fivebar1TE_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin);
Ja_rot = fivebar1TE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
end