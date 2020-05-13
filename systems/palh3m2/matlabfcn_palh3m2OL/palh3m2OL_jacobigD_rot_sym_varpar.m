% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m2OL
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% JgD_rot [3x10]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh3m2OL_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),uint8(0),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_jacobigD_rot_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_jacobigD_rot_sym_varpar: qJD has to be [10x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2OL_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_jacobigD_rot_sym_varpar: pkin has to be [16x1] (double)');
JgD_rot=NaN(3,10);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0, 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = qJD(1) * cos(qJ(1));
	t35 = qJD(1) * sin(qJ(1));
	t1 = [0, t36, t36, 0, 0, 0, 0, 0, 0, 0; 0, t35, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (6->2), div. (0->0), fcn. (6->2), ass. (0->3)
	t41 = qJD(1) * cos(qJ(1));
	t40 = qJD(1) * sin(qJ(1));
	t1 = [0, t41, t41, t41, 0, 0, 0, 0, 0, 0; 0, t40, t40, t40, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:52
	% EndTime: 2020-05-07 04:44:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->5), mult. (15->8), div. (0->0), fcn. (15->4), ass. (0->9)
	t114 = qJD(2) + qJD(3) + qJD(4);
	t117 = qJ(2) + qJ(3) + qJ(4);
	t120 = cos(t117) * t114;
	t118 = sin(qJ(1));
	t115 = qJD(1) * t118;
	t119 = cos(qJ(1));
	t116 = qJD(1) * t119;
	t112 = sin(t117);
	t1 = [0, t116, t116, t116, t112 * t115 - t119 * t120, 0, 0, 0, 0, 0; 0, t115, t115, t115, -t112 * t116 - t118 * t120, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t114 * t112, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t28 = qJD(1) * cos(qJ(1));
	t27 = qJD(1) * sin(qJ(1));
	t1 = [0, t28, 0, 0, 0, 0, t28, 0, 0, 0; 0, t27, 0, 0, 0, 0, t27, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobigD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (6->2), div. (0->0), fcn. (6->2), ass. (0->3)
	t35 = qJD(1) * cos(qJ(1));
	t34 = qJD(1) * sin(qJ(1));
	t1 = [0, t35, 0, 0, 0, 0, t35, t35, 0, 0; 0, t34, 0, 0, 0, 0, t34, t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
end