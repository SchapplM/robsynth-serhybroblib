% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh4m1OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% JR_rot [9x8]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh4m1OL_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_jacobiR_rot_sym_varpar: qJ has to be [8x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1OL_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->72)
	unknown=NaN(9,8);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(4,1) = 0;
	unknown(4,2) = 0;
	unknown(4,3) = 0;
	unknown(4,4) = 0;
	unknown(4,5) = 0;
	unknown(4,6) = 0;
	unknown(4,7) = 0;
	unknown(4,8) = 0;
	unknown(5,1) = 0;
	unknown(5,2) = 0;
	unknown(5,3) = 0;
	unknown(5,4) = 0;
	unknown(5,5) = 0;
	unknown(5,6) = 0;
	unknown(5,7) = 0;
	unknown(5,8) = 0;
	unknown(6,1) = 0;
	unknown(6,2) = 0;
	unknown(6,3) = 0;
	unknown(6,4) = 0;
	unknown(6,5) = 0;
	unknown(6,6) = 0;
	unknown(6,7) = 0;
	unknown(6,8) = 0;
	unknown(7,1) = 0;
	unknown(7,2) = 0;
	unknown(7,3) = 0;
	unknown(7,4) = 0;
	unknown(7,5) = 0;
	unknown(7,6) = 0;
	unknown(7,7) = 0;
	unknown(7,8) = 0;
	unknown(8,1) = 0;
	unknown(8,2) = 0;
	unknown(8,3) = 0;
	unknown(8,4) = 0;
	unknown(8,5) = 0;
	unknown(8,6) = 0;
	unknown(8,7) = 0;
	unknown(8,8) = 0;
	unknown(9,1) = 0;
	unknown(9,2) = 0;
	unknown(9,3) = 0;
	unknown(9,4) = 0;
	unknown(9,5) = 0;
	unknown(9,6) = 0;
	unknown(9,7) = 0;
	unknown(9,8) = 0;
	JR_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->74)
	unknown=NaN(9,8);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = -t1;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t2;
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t2;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = -t1;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = 0.0e0;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = 0.0e0;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->80)
	unknown=NaN(9,8);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = t1 * t2;
	t4 = cos(qJ(1));
	t5 = sin(qJ(2));
	t6 = t4 * t5;
	t7 = t4 * t2;
	t8 = t1 * t5;
	unknown(1,1) = -t3;
	unknown(1,2) = -t6;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t7;
	unknown(2,2) = -t8;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t2;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = t8;
	unknown(4,2) = -t7;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = -t6;
	unknown(5,2) = -t3;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t5;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = t4;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->80)
	unknown=NaN(9,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = cos(qJ(1));
	t5 = cos(qJ(2));
	t6 = t4 * t5;
	t7 = t4 * t2;
	t8 = t1 * t5;
	unknown(1,1) = -t3;
	unknown(1,2) = t6;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t7;
	unknown(2,2) = t8;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t2;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t4;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = -t1;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = -t8;
	unknown(7,2) = -t7;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t6;
	unknown(8,2) = -t3;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t5;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (16->9), mult. (56->16), div. (0->0), fcn. (90->6), ass. (0->88)
	unknown=NaN(9,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = sin(qJ(4));
	t6 = cos(qJ(2));
	t7 = t1 * t6;
	t8 = cos(qJ(4));
	t10 = t3 * t4 - t7 * t8;
	t11 = cos(qJ(1));
	t12 = t11 * t2;
	t14 = t11 * t6;
	t16 = -t12 * t8 - t14 * t4;
	t19 = -t12 * t4 + t14 * t8;
	t22 = -t3 * t8 - t7 * t4;
	t25 = -t2 * t4 + t6 * t8;
	t28 = -t2 * t8 - t6 * t4;
	unknown(1,1) = t10;
	unknown(1,2) = t16;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = t16;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t19;
	unknown(2,2) = t22;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = t22;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t25;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = t25;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t22;
	unknown(4,2) = -t19;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = -t19;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t16;
	unknown(5,2) = t10;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = t10;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t28;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = t28;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = t11;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (66->11), mult. (76->16), div. (0->0), fcn. (122->6), ass. (0->89)
	unknown=NaN(9,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = qJ(4) + qJ(5);
	t5 = sin(t4);
	t7 = cos(qJ(2));
	t8 = t1 * t7;
	t9 = cos(t4);
	t11 = t3 * t5 - t8 * t9;
	t12 = cos(qJ(1));
	t13 = t12 * t2;
	t15 = t12 * t7;
	t17 = -t13 * t9 - t15 * t5;
	t20 = -t13 * t5 + t15 * t9;
	t23 = -t3 * t9 - t8 * t5;
	t26 = -t2 * t5 + t7 * t9;
	t29 = -t2 * t9 - t7 * t5;
	unknown(1,1) = t11;
	unknown(1,2) = t17;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = t17;
	unknown(1,5) = t17;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t20;
	unknown(2,2) = t23;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = t23;
	unknown(2,5) = t23;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t26;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = t26;
	unknown(3,5) = t26;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t23;
	unknown(4,2) = -t20;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = -t20;
	unknown(4,5) = -t20;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t17;
	unknown(5,2) = t11;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = t11;
	unknown(5,5) = t11;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t29;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = t29;
	unknown(6,5) = t29;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = t12;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:11
	% EndTime: 2020-04-11 23:04:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (134->27), mult. (170->34), div. (0->0), fcn. (256->8), ass. (0->101)
	unknown=NaN(9,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = qJ(4) + qJ(5);
	t5 = cos(t4);
	t7 = cos(qJ(2));
	t8 = t1 * t7;
	t9 = sin(t4);
	t11 = t3 * t5 + t8 * t9;
	t12 = cos(qJ(6));
	t14 = cos(qJ(1));
	t15 = sin(qJ(6));
	t16 = t14 * t15;
	t18 = t14 * t7;
	t20 = t14 * t2;
	t22 = -t18 * t5 + t20 * t9;
	t23 = t22 * t12;
	t26 = -t18 * t9 - t20 * t5;
	t29 = -t1 * t12 + t26 * t15;
	t32 = -t1 * t15 - t26 * t12;
	t35 = t3 * t9 - t8 * t5;
	t36 = t35 * t12;
	t38 = t14 * t12;
	t42 = -t2 * t5 - t7 * t9;
	t43 = t42 * t12;
	t46 = -t2 * t9 + t7 * t5;
	t50 = t22 * t15;
	t51 = t35 * t15;
	t54 = t42 * t15;
	unknown(1,1) = -t11 * t12 - t16;
	unknown(1,2) = -t23;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = -t23;
	unknown(1,5) = -t23;
	unknown(1,6) = t29;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t32;
	unknown(2,2) = -t36;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = -t36;
	unknown(2,5) = -t36;
	unknown(2,6) = -t11 * t15 + t38;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t43;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = -t43;
	unknown(3,5) = -t43;
	unknown(3,6) = t46 * t15;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = t11 * t15 - t38;
	unknown(4,2) = t50;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = t50;
	unknown(4,5) = t50;
	unknown(4,6) = -t32;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t29;
	unknown(5,2) = t51;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = t51;
	unknown(5,5) = t51;
	unknown(5,6) = -t11 * t12 - t16;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t54;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = t54;
	unknown(6,5) = t54;
	unknown(6,6) = t46 * t12;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = t35;
	unknown(7,2) = t26;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = t26;
	unknown(7,5) = t26;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = -t22;
	unknown(8,2) = -t11;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = -t11;
	unknown(8,5) = -t11;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t46;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = t46;
	unknown(9,5) = t46;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:12
	% EndTime: 2020-04-11 23:04:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->80)
	unknown=NaN(9,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(7));
	t3 = t1 * t2;
	t4 = cos(qJ(1));
	t5 = cos(qJ(7));
	t6 = t4 * t5;
	t7 = t4 * t2;
	t8 = t1 * t5;
	unknown(1,1) = t3;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = -t6;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t7;
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = -t8;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = -t2;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = t8;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = t7;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = -t6;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = t3;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = -t5;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = t4;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JR_rot = unknown;
else
	JR_rot=NaN(9,8);
end