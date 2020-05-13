% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% JR_rot [9x12]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = picker2Dm2OL_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_jacobiR_rot_sym_varpar: qJ has to be [12x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'picker2Dm2OL_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
JR_rot=NaN(9,12);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->3), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t18 = qJ(1) + qJ(2);
	t17 = cos(t18);
	t16 = sin(t18);
	t1 = [t16, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t17, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t17, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t16, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t21 = qJ(1) + qJ(2) + qJ(3);
	t20 = cos(t21);
	t19 = sin(t21);
	t1 = [-t19, -t19, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0; t20, t20, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t20, -t20, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t19, -t19, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (27->4), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t23 = qJ(1) + qJ(2) + qJ(4);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21, t21, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0; -t22, -t22, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t22, t22, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0; t21, t21, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = pkin(8) + qJ(5);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t16 = qJ(1) + qJ(2) + qJ(6);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0; t15, t15, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t15, -t15, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0; -t14, -t14, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(7));
	t8 = sin(qJ(7));
	t1 = [0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(8);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t15, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiR_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (52->5), mult. (0->0), div. (0->0), fcn. (16->2), ass. (0->4)
	t28 = qJ(1) + qJ(2) + qJ(3) + qJ(9);
	t27 = cos(t28);
	t26 = sin(t28);
	t1 = [t26, t26, t26, 0, 0, 0, 0, 0, t26, 0, 0, 0; -t27, -t27, -t27, 0, 0, 0, 0, 0, -t27, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t27, t27, t27, 0, 0, 0, 0, 0, t27, 0, 0, 0; t26, t26, t26, 0, 0, 0, 0, 0, t26, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiR_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:38
	% EndTime: 2020-05-09 23:20:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (60->13), mult. (0->0), div. (0->0), fcn. (16->2), ass. (0->4)
	t26 = qJ(1) + qJ(2) + qJ(4) + qJ(10);
	t25 = cos(t26);
	t24 = sin(t26);
	t1 = [-t24, -t24, 0, -t24, 0, 0, 0, 0, 0, -t24, 0, 0; t25, t25, 0, t25, 0, 0, 0, 0, 0, t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t25, -t25, 0, -t25, 0, 0, 0, 0, 0, -t25, 0, 0; -t24, -t24, 0, -t24, 0, 0, 0, 0, 0, -t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
end