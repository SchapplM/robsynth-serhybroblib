% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% JR_rot [9x1]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = fourbar1TE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_jacobiR_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1TE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_jacobiR_rot_sym_varpar: pkin has to be [4x1] (double)');
JR_rot=NaN(9,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8; t9; 0; -t9; -t8; 0; 0; 0; 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:43
	% EndTime: 2020-04-24 19:49:43
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (284->26), mult. (400->40), div. (16->4), fcn. (110->4), ass. (0->26)
	t93 = sin(qJ(1));
	t107 = pkin(2) * t93;
	t94 = cos(qJ(1));
	t106 = pkin(2) * t94;
	t104 = (-0.2e1 * t106 + pkin(1)) * pkin(1);
	t96 = pkin(2) ^ 2;
	t89 = t96 + t104;
	t86 = pkin(3) ^ 2 - pkin(4) ^ 2 + t89;
	t83 = t86 * t107;
	t87 = 0.1e1 / t89;
	t88 = 0.1e1 / t89 ^ 2;
	t90 = -pkin(1) + t106;
	t109 = -pkin(3) - pkin(4);
	t84 = (pkin(2) - t109) * (pkin(2) + t109) + t104;
	t108 = -pkin(3) + pkin(4);
	t85 = (pkin(2) - t108) * (pkin(2) + t108) + t104;
	t98 = sqrt(-t84 * t85);
	t112 = ((t90 * t98 + t83) * t88 * t107 - t96 * t93 ^ 2 * t87) * pkin(1);
	t110 = t87 / 0.2e1;
	t103 = pkin(1) * t107;
	t105 = (-t84 - t85) * t103 / t98;
	t80 = t98 * t107;
	t99 = t90 * t105 + t86 * t106 - t80;
	t95 = 0.1e1 / pkin(3);
	t78 = ((t83 + (t93 * t105 + t94 * t98) * pkin(2)) * t110 + (-t90 * t87 - (-t90 * t86 + t80) * t88) * t103) * t95;
	t1 = [t78; (-t99 * t87 / 0.2e1 + t112) * t95; 0; (t99 * t110 - t112) * t95; t78; 0; 0; 0; 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:43
	% EndTime: 2020-04-24 19:49:43
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (284->26), mult. (400->41), div. (16->4), fcn. (110->4), ass. (0->25)
	t92 = sin(qJ(1));
	t106 = pkin(2) * t92;
	t93 = cos(qJ(1));
	t105 = pkin(2) * t93;
	t102 = (-0.2e1 * t105 + pkin(1)) * pkin(1);
	t95 = pkin(2) ^ 2;
	t88 = t95 + t102;
	t85 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t88;
	t86 = 0.1e1 / t88;
	t87 = 0.1e1 / t88 ^ 2;
	t89 = -pkin(1) + t105;
	t108 = -pkin(3) - pkin(4);
	t83 = (pkin(2) - t108) * (pkin(2) + t108) + t102;
	t107 = -pkin(3) + pkin(4);
	t84 = (pkin(2) - t107) * (pkin(2) + t107) + t102;
	t97 = sqrt(-t83 * t84);
	t111 = ((t85 * t106 - t89 * t97) * t87 * t106 - t95 * t92 ^ 2 * t86) * pkin(1);
	t109 = t86 / 0.2e1;
	t104 = t92 * pkin(1);
	t103 = 0.1e1 / t97 * (-t83 - t84) * pkin(2) * t104;
	t80 = t97 * t106;
	t98 = -t89 * t103 + t85 * t105 + t80;
	t94 = 0.1e1 / pkin(4);
	t78 = ((t93 * t97 + (-t85 + t103) * t92) * t109 + (t89 * t86 - (t89 * t85 + t80) * t87) * t104) * t94 * pkin(2);
	t1 = [t78; (t98 * t109 - t111) * t94; 0; (-t98 * t86 / 0.2e1 + t111) * t94; t78; 0; 0; 0; 0;];
	JR_rot = t1;
end