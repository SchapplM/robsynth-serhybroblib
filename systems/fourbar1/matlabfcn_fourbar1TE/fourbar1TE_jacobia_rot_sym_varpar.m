% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in fourbar1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% Ja_rot [3x1]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = fourbar1TE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_jacobia_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1TE_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_jacobia_rot_sym_varpar: pkin has to be [4x1] (double)');
Ja_rot=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0; 0; 1;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:43
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (282->27), mult. (383->45), div. (14->6), fcn. (106->4), ass. (0->24)
	t65 = -pkin(3) - pkin(4);
	t64 = -pkin(3) + pkin(4);
	t51 = sin(qJ(1));
	t63 = pkin(2) * t51;
	t52 = cos(qJ(1));
	t62 = pkin(2) * t52;
	t59 = (-0.2e1 * t62 + pkin(1)) * pkin(1);
	t43 = (pkin(2) - t65) * (pkin(2) + t65) + t59;
	t44 = (pkin(2) - t64) * (pkin(2) + t64) + t59;
	t56 = sqrt(-t43 * t44);
	t58 = pkin(1) * t63;
	t61 = (-t43 - t44) * t58 / t56;
	t60 = t51 * t56;
	t54 = pkin(2) ^ 2;
	t48 = t54 + t59;
	t49 = -pkin(1) + t62;
	t47 = 0.1e1 / t48 ^ 2;
	t46 = 0.1e1 / t48;
	t45 = pkin(3) ^ 2 - pkin(4) ^ 2 + t48;
	t42 = t45 * t63;
	t38 = t49 * t56 + t42;
	t37 = pkin(2) * t60 - t45 * t49;
	t36 = 0.1e1 / t37 ^ 2;
	t1 = [0; 0; 0.2e1 * ((-(t49 * t61 + (t45 * t52 - t60) * pkin(2)) * t46 / 0.2e1 + (pkin(2) * t38 * t47 - t46 * t51 * t54) * t51 * pkin(1)) / t37 + ((t42 + (t51 * t61 + t52 * t56) * pkin(2)) * t46 / 0.2e1 + (-t37 * t47 - t46 * t49) * t58) * t38 * t36) / (t36 * t38 ^ 2 + 0.1e1) * t48;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:43
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (282->27), mult. (383->44), div. (14->6), fcn. (106->4), ass. (0->24)
	t65 = -pkin(3) - pkin(4);
	t64 = -pkin(3) + pkin(4);
	t50 = sin(qJ(1));
	t63 = pkin(1) * t50;
	t53 = pkin(2) ^ 2;
	t51 = cos(qJ(1));
	t60 = pkin(2) * t51;
	t57 = (-0.2e1 * t60 + pkin(1)) * pkin(1);
	t47 = t53 + t57;
	t44 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t47;
	t48 = -pkin(1) + t60;
	t42 = (pkin(2) - t65) * (pkin(2) + t65) + t57;
	t43 = (pkin(2) - t64) * (pkin(2) + t64) + t57;
	t55 = sqrt(-t42 * t43);
	t61 = pkin(2) * t50;
	t38 = -t44 * t61 + t48 * t55;
	t62 = pkin(2) * t38;
	t59 = 0.1e1 / t55 * (-t42 - t43) * pkin(1) * t61;
	t58 = t50 * t55;
	t46 = 0.1e1 / t47 ^ 2;
	t45 = 0.1e1 / t47;
	t37 = pkin(2) * t58 + t44 * t48;
	t36 = 0.1e1 / t37 ^ 2;
	t1 = [0; 0; 0.2e1 * ((-(t48 * t59 + (-t51 * t44 - t58) * pkin(2)) * t45 / 0.2e1 + (t45 * t50 * t53 + t46 * t62) * t63) / t37 + ((t51 * t55 + (-t44 + t59) * t50) * t45 / 0.2e1 + (-t37 * t46 + t45 * t48) * t63) * t36 * t62) / (t36 * t38 ^ 2 + 0.1e1) * t47;];
	Ja_rot = t1;
end