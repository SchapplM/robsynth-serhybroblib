% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JgD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh3m2DE1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_jacobigD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_jacobigD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_jacobigD_rot_sym_varpar: pkin has to be [18x1] (double)');
JgD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:34
	% EndTime: 2020-05-07 01:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:34
	% EndTime: 2020-05-07 01:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:35
	% EndTime: 2020-05-07 01:57:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:35
	% EndTime: 2020-05-07 01:57:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = qJD(1) * cos(qJ(1));
	t42 = qJD(1) * sin(qJ(1));
	t1 = [0, t43, t43, 0; 0, t42, t42, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:36
	% EndTime: 2020-05-07 01:57:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:36
	% EndTime: 2020-05-07 01:57:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->5), mult. (16->9), div. (0->0), fcn. (22->8), ass. (0->7)
	t116 = pkin(17) + pkin(18);
	t117 = sin(pkin(16));
	t118 = cos(pkin(16));
	t119 = sin(pkin(15));
	t120 = cos(pkin(15));
	t121 = qJD(1) * ((t117 * t120 + t118 * t119) * cos(t116) - sin(t116) * (t117 * t119 - t118 * t120));
	t1 = [0, 0, 0, -sin(qJ(1)) * t121; 0, 0, 0, cos(qJ(1)) * t121; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:37
	% EndTime: 2020-05-07 01:57:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:38
	% EndTime: 2020-05-07 01:57:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobigD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:39
	% EndTime: 2020-05-07 01:57:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t41 = qJD(1) * cos(qJ(1));
	t40 = qJD(1) * sin(qJ(1));
	t1 = [0, t41, t41, 0; 0, t40, t40, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
end