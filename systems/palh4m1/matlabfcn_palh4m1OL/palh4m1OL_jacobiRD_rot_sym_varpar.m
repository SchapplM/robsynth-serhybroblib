% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh4m1OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% qJD [8x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% JRD_rot [9x8]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh4m1OL_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(8,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_jacobiRD_rot_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [8 1]), ...
  'palh4m1OL_jacobiRD_rot_sym_varpar: qJD has to be [8x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1OL_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_jacobiRD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
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
	JRD_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->76)
	unknown=NaN(9,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	unknown(1,1) = -t2;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t4;
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
	unknown(4,1) = t4;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = -t2;
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
	JRD_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->84)
	unknown=NaN(9,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = cos(qJ(2));
	t5 = sin(qJ(1));
	t6 = t5 * qJD(2);
	t7 = sin(qJ(2));
	t9 = -t2 * t3 + t6 * t7;
	t10 = qJD(1) * t5;
	t12 = t1 * qJD(2);
	t14 = t10 * t7 - t12 * t3;
	t17 = -t10 * t3 - t12 * t7;
	t20 = -t2 * t7 - t6 * t3;
	unknown(1,1) = t9;
	unknown(1,2) = t14;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t17;
	unknown(2,2) = t20;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -qJD(2) * t7;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t20;
	unknown(4,2) = -t17;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t14;
	unknown(5,2) = t9;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -qJD(2) * t3;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = -t10;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t2;
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
	JRD_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->84)
	unknown=NaN(9,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(2));
	t5 = sin(qJ(1));
	t6 = t5 * qJD(2);
	t7 = cos(qJ(2));
	t9 = -t2 * t3 - t6 * t7;
	t10 = qJD(1) * t5;
	t12 = t1 * qJD(2);
	t14 = -t10 * t7 - t12 * t3;
	t17 = -t10 * t3 + t12 * t7;
	t20 = t2 * t7 - t6 * t3;
	unknown(1,1) = t9;
	unknown(1,2) = t14;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t17;
	unknown(2,2) = t20;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = qJD(2) * t7;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = t10;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = -t2;
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
	unknown(7,1) = -t20;
	unknown(7,2) = -t17;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t14;
	unknown(8,2) = t9;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = -qJD(2) * t3;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JRD_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (250->50), div. (0->0), fcn. (250->6), ass. (0->102)
	unknown=NaN(9,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(2));
	t4 = sin(qJ(4));
	t5 = t3 * t4;
	t7 = sin(qJ(1));
	t8 = t7 * qJD(2);
	t9 = cos(qJ(2));
	t10 = t9 * t4;
	t12 = t7 * t3;
	t13 = cos(qJ(4));
	t14 = qJD(4) * t13;
	t16 = t9 * t13;
	t18 = t3 * t13;
	t20 = t7 * t9;
	t21 = qJD(4) * t4;
	t23 = t8 * t10 + t12 * t14 - t2 * t16 + t8 * t18 + t2 * t5 + t20 * t21;
	t24 = qJD(1) * t7;
	t26 = t1 * qJD(2);
	t28 = t1 * t3;
	t32 = t1 * t9;
	t34 = t24 * t10 - t32 * t14 - t26 * t16 + t24 * t18 + t28 * t21 + t26 * t5;
	t41 = -t26 * t10 - t28 * t14 - t24 * t16 - t26 * t18 - t32 * t21 + t24 * t5;
	t48 = -t2 * t10 + t12 * t21 - t20 * t14 - t8 * t16 - t2 * t18 + t8 * t5;
	t49 = qJD(2) * t3;
	t51 = t9 * qJD(4);
	t53 = qJD(2) * t9;
	t55 = t3 * qJD(4);
	t57 = -t49 * t13 - t55 * t13 - t51 * t4 - t53 * t4;
	t62 = -t51 * t13 - t53 * t13 + t49 * t4 + t55 * t4;
	unknown(1,1) = t23;
	unknown(1,2) = t34;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = t34;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t41;
	unknown(2,2) = t48;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = t48;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t57;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = t57;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t48;
	unknown(4,2) = -t41;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = -t41;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t34;
	unknown(5,2) = t23;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = t23;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t62;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = t62;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = -t24;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t2;
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
	JRD_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:11
	% EndTime: 2020-04-11 23:04:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (263->33), mult. (338->50), div. (0->0), fcn. (338->6), ass. (0->104)
	unknown=NaN(9,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(2));
	t4 = qJ(4) + qJ(5);
	t5 = sin(t4);
	t6 = t3 * t5;
	t8 = sin(qJ(1));
	t9 = t8 * qJD(2);
	t10 = cos(qJ(2));
	t11 = t10 * t5;
	t13 = t8 * t3;
	t14 = qJD(4) + qJD(5);
	t15 = cos(t4);
	t16 = t14 * t15;
	t18 = t10 * t15;
	t20 = t3 * t15;
	t22 = t8 * t10;
	t23 = t14 * t5;
	t25 = t9 * t11 + t13 * t16 - t2 * t18 + t2 * t6 + t9 * t20 + t22 * t23;
	t26 = qJD(1) * t8;
	t28 = t1 * qJD(2);
	t30 = t1 * t3;
	t34 = t1 * t10;
	t36 = t26 * t11 - t34 * t16 - t28 * t18 + t26 * t20 + t30 * t23 + t28 * t6;
	t43 = -t28 * t11 - t30 * t16 - t26 * t18 - t28 * t20 - t34 * t23 + t26 * t6;
	t50 = -t2 * t11 + t13 * t23 - t22 * t16 - t9 * t18 - t2 * t20 + t9 * t6;
	t51 = qJD(2) * t3;
	t53 = t10 * t14;
	t55 = qJD(2) * t10;
	t57 = t3 * t14;
	t59 = -t51 * t15 - t57 * t15 - t53 * t5 - t55 * t5;
	t64 = -t53 * t15 - t55 * t15 + t51 * t5 + t57 * t5;
	unknown(1,1) = t25;
	unknown(1,2) = t36;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = t36;
	unknown(1,5) = t36;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t43;
	unknown(2,2) = t50;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = t50;
	unknown(2,5) = t50;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t59;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = t59;
	unknown(3,5) = t59;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t50;
	unknown(4,2) = -t43;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = -t43;
	unknown(4,5) = -t43;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t36;
	unknown(5,2) = t25;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = t25;
	unknown(5,5) = t25;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t64;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = t64;
	unknown(6,5) = t64;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = -t26;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t2;
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
	JRD_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:11
	% EndTime: 2020-04-11 23:04:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (591->66), mult. (804->103), div. (0->0), fcn. (832->8), ass. (0->128)
	unknown=NaN(9,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(2));
	t4 = qJ(4) + qJ(5);
	t5 = cos(t4);
	t6 = t3 * t5;
	t8 = sin(qJ(1));
	t9 = t8 * qJD(2);
	t10 = cos(qJ(2));
	t11 = t10 * t5;
	t13 = t8 * t3;
	t14 = qJD(4) + qJD(5);
	t15 = sin(t4);
	t16 = t14 * t15;
	t18 = t10 * t15;
	t20 = t3 * t15;
	t22 = t8 * t10;
	t23 = t14 * t5;
	t25 = t9 * t11 - t13 * t16 + t2 * t18 + t2 * t6 - t9 * t20 + t22 * t23;
	t26 = cos(qJ(6));
	t30 = t13 * t5 + t22 * t15;
	t31 = t30 * qJD(6);
	t32 = sin(qJ(6));
	t34 = qJD(1) * t8;
	t35 = t34 * t32;
	t36 = t1 * qJD(6);
	t37 = t36 * t26;
	t40 = t1 * qJD(2);
	t42 = t1 * t10;
	t46 = t1 * t3;
	t48 = t34 * t11 + t42 * t16 + t40 * t18 - t34 * t20 + t46 * t23 + t40 * t6;
	t53 = (t46 * t15 - t42 * t5) * qJD(6);
	t55 = -t48 * t26 + t53 * t32;
	t62 = -t40 * t11 + t46 * t16 + t34 * t18 + t40 * t20 - t42 * t23 + t34 * t6;
	t67 = (-t42 * t15 - t46 * t5) * qJD(6);
	t70 = t8 * qJD(6);
	t72 = -t2 * t26 + t67 * t26 + t62 * t32 + t70 * t32;
	t77 = -t2 * t32 - t62 * t26 - t70 * t26 + t67 * t32;
	t84 = -t2 * t11 + t13 * t23 + t22 * t16 + t9 * t18 + t2 * t20 + t9 * t6;
	t89 = (t13 * t15 - t22 * t5) * qJD(6);
	t91 = -t84 * t26 + t89 * t32;
	t93 = -t30 * qJD(6);
	t95 = t34 * t26;
	t96 = t36 * t32;
	t98 = qJD(2) * t10;
	t100 = t3 * t14;
	t102 = qJD(2) * t3;
	t104 = t10 * t14;
	t106 = t100 * t15 + t102 * t15 - t104 * t5 - t98 * t5;
	t109 = (-t6 - t18) * qJD(6);
	t111 = -t106 * t26 + t109 * t32;
	t116 = -t100 * t5 - t102 * t5 - t104 * t15 - t98 * t15;
	t119 = (t11 - t20) * qJD(6);
	t127 = t53 * t26 + t48 * t32;
	t130 = t89 * t26 + t84 * t32;
	t136 = t106 * t32 + t109 * t26;
	unknown(1,1) = -t25 * t26 + t31 * t32 + t35 - t37;
	unknown(1,2) = t55;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = t55;
	unknown(1,5) = t55;
	unknown(1,6) = t72;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t77;
	unknown(2,2) = t91;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = t91;
	unknown(2,5) = t91;
	unknown(2,6) = -t25 * t32 + t93 * t26 - t95 - t96;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t111;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = t111;
	unknown(3,5) = t111;
	unknown(3,6) = t116 * t32 + t119 * t26;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = t25 * t32 + t31 * t26 + t95 + t96;
	unknown(4,2) = t127;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = t127;
	unknown(4,5) = t127;
	unknown(4,6) = -t77;
	unknown(4,7) = 0.0e0;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t72;
	unknown(5,2) = t130;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = t130;
	unknown(5,5) = t130;
	unknown(5,6) = -t25 * t26 - t93 * t32 + t35 - t37;
	unknown(5,7) = 0.0e0;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t136;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = t136;
	unknown(6,5) = t136;
	unknown(6,6) = t116 * t26 - t119 * t32;
	unknown(6,7) = 0.0e0;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = t84;
	unknown(7,2) = t62;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = t62;
	unknown(7,5) = t62;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = -t48;
	unknown(8,2) = -t25;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = -t25;
	unknown(8,5) = -t25;
	unknown(8,6) = 0.0e0;
	unknown(8,7) = 0.0e0;
	unknown(8,8) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t116;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = t116;
	unknown(9,5) = t116;
	unknown(9,6) = 0.0e0;
	unknown(9,7) = 0.0e0;
	unknown(9,8) = 0.0e0;
	JRD_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:12
	% EndTime: 2020-04-11 23:04:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->84)
	unknown=NaN(9,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(7));
	t5 = sin(qJ(1));
	t6 = t5 * qJD(7);
	t7 = cos(qJ(7));
	t9 = t2 * t3 + t6 * t7;
	t10 = qJD(1) * t5;
	t12 = t1 * qJD(7);
	t14 = t10 * t7 + t12 * t3;
	t17 = t10 * t3 - t12 * t7;
	t20 = -t2 * t7 + t6 * t3;
	unknown(1,1) = t9;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = t14;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t17;
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = t20;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = -qJD(7) * t7;
	unknown(3,8) = 0.0e0;
	unknown(4,1) = -t20;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(4,6) = 0.0e0;
	unknown(4,7) = -t17;
	unknown(4,8) = 0.0e0;
	unknown(5,1) = t14;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(5,6) = 0.0e0;
	unknown(5,7) = t9;
	unknown(5,8) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(6,6) = 0.0e0;
	unknown(6,7) = qJD(7) * t3;
	unknown(6,8) = 0.0e0;
	unknown(7,1) = -t10;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(7,6) = 0.0e0;
	unknown(7,7) = 0.0e0;
	unknown(7,8) = 0.0e0;
	unknown(8,1) = t2;
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
	JRD_rot = unknown;
else
	JRD_rot=NaN(9,8);
end