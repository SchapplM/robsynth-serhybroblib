% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fourbar1TE
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
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
% Jg_rot [3x1]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = fourbar1TE_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_jacobig_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1TE_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_jacobig_rot_sym_varpar: pkin has to be [4x1] (double)');
Jg_rot=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 1;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:43
	% EndTime: 2020-04-24 19:49:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->16), mult. (42->14), div. (1->2), fcn. (12->4), ass. (0->7)
	t77 = -pkin(3) - pkin(4);
	t76 = -pkin(3) + pkin(4);
	t75 = pkin(2) * cos(qJ(1));
	t74 = (-0.2e1 * t75 + pkin(1)) * pkin(1);
	t73 = pkin(2) ^ 2 + t74;
	t72 = sqrt(-((pkin(2) - t77) * (pkin(2) + t77) + t74) * ((pkin(2) - t76) * (pkin(2) + t76) + t74));
	t1 = [0; 0; 0.1e1 - ((pkin(1) - t75) * t72 + pkin(2) * sin(qJ(1)) * (-pkin(3) ^ 2 + pkin(4) ^ 2 + t73)) / t72 * pkin(1) / t73;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:43
	% EndTime: 2020-04-24 19:49:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (35->16), mult. (42->15), div. (1->2), fcn. (12->4), ass. (0->7)
	t77 = -pkin(3) - pkin(4);
	t76 = -pkin(3) + pkin(4);
	t75 = pkin(1) * cos(qJ(1));
	t74 = pkin(1) ^ 2 - 0.2e1 * pkin(2) * t75;
	t73 = pkin(2) ^ 2 + t74;
	t72 = sqrt(-((pkin(2) - t77) * (pkin(2) + t77) + t74) * ((pkin(2) - t76) * (pkin(2) + t76) + t74));
	t1 = [0; 0; -((-pkin(2) + t75) * t72 - pkin(1) * sin(qJ(1)) * (pkin(3) ^ 2 - pkin(4) ^ 2 + t73)) * pkin(2) / t72 / t73;];
	Jg_rot = t1;
end