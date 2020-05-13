% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh4m1OL
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
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
% JgD_rot [3x8]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh4m1OL_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(8,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_jacobigD_rot_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [8 1]), ...
  'palh4m1OL_jacobigD_rot_sym_varpar: qJD has to be [8x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1OL_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_jacobigD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->24)
	unknown=NaN(3,8);
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
	JgD_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->24)
	unknown=NaN(3,8);
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
	JgD_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->26)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t3 = sin(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (qJD(1) * t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = (qJD(1) * t3);
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
	JgD_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->26)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t3 = sin(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (qJD(1) * t1);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = (qJD(1) * t3);
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
	JgD_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->28)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	unknown(1,1) = 0;
	unknown(1,2) = t2;
	unknown(1,3) = 0;
	unknown(1,4) = t2;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = t4;
	unknown(2,3) = 0;
	unknown(2,4) = t4;
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
	JgD_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:11
	% EndTime: 2020-04-11 23:04:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (6->2), div. (0->0), fcn. (6->2), ass. (0->28)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	unknown(1,1) = 0;
	unknown(1,2) = t2;
	unknown(1,3) = 0;
	unknown(1,4) = t2;
	unknown(1,5) = t2;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = t4;
	unknown(2,3) = 0;
	unknown(2,4) = t4;
	unknown(2,5) = t4;
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
	JgD_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:11
	% EndTime: 2020-04-11 23:04:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (35->15), mult. (50->34), div. (0->0), fcn. (50->6), ass. (0->42)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t2 = qJD(1) * t1;
	t3 = sin(qJ(1));
	t4 = qJD(1) * t3;
	t5 = sin(qJ(2));
	t6 = qJ(4) + qJ(5);
	t7 = sin(t6);
	t8 = t7 * t5;
	t10 = t1 * qJD(2);
	t11 = cos(qJ(2));
	t12 = t11 * t7;
	t15 = qJD(4) + qJD(5);
	t16 = cos(t6);
	t17 = t15 * t16;
	t19 = t11 * t16;
	t21 = t5 * t16;
	t24 = t15 * t7;
	t28 = t3 * qJD(2);
	unknown(1,1) = 0;
	unknown(1,2) = t2;
	unknown(1,3) = 0;
	unknown(1,4) = t2;
	unknown(1,5) = t2;
	unknown(1,6) = (-t1 * t11 * t24 - t1 * t5 * t17 - t10 * t12 - t10 * t21 - t4 * t19 + t4 * t8);
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = t4;
	unknown(2,3) = 0;
	unknown(2,4) = t4;
	unknown(2,5) = t4;
	unknown(2,6) = (-t3 * t11 * t24 - t3 * t5 * t17 - t28 * t12 + t2 * t19 - t2 * t8 - t28 * t21);
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = (qJD(2) * t11 * t16 - qJD(2) * t5 * t7 + t11 * t15 * t16 - t5 * t15 * t7);
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	JgD_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:12
	% EndTime: 2020-04-11 23:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->26)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t3 = sin(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = (qJD(1) * t1);
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = (qJD(1) * t3);
	unknown(2,8) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	JgD_rot = unknown;
else
	JgD_rot=NaN(3,8);
end