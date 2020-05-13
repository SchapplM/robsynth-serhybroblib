% Analytische Jacobi-Matrix für beliebiges Segment von
% mg10hlOL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in mg10hlOL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% Ja [6x13]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:06
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = mg10hlOL_jacobia_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),uint8(0),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlOL_jacobia_sym_varpar: qJ has to be [13x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'mg10hlOL_jacobia_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlOL_jacobia_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlOL_jacobia_sym_varpar: pkin has to be [16x1] (double)');

% Function calls
Ja_transl = mg10hlOL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin);
Ja_rot = mg10hlOL_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
end