% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m2DE2
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
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh3m2DE2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_jacobigD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_jacobigD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_jacobigD_rot_sym_varpar: pkin has to be [18x1] (double)');
JgD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:22
	% EndTime: 2020-05-07 02:13:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = qJD(1) * cos(qJ(1));
	t35 = qJD(1) * sin(qJ(1));
	t1 = [0, t36, t36, 0; 0, t35, t35, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:25
	% EndTime: 2020-05-07 02:13:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:32
	% EndTime: 2020-05-07 02:13:32
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (753->24), mult. (1377->40), div. (18->3), fcn. (1980->15), ass. (0->32)
	t330 = sin(pkin(16));
	t331 = cos(pkin(16));
	t335 = sin(pkin(15));
	t339 = cos(pkin(15));
	t325 = t330 * t339 + t331 * t335;
	t326 = -t330 * t335 + t331 * t339;
	t332 = sin(qJ(3));
	t336 = cos(qJ(3));
	t340 = t325 * t332 - t326 * t336;
	t320 = t340 * qJD(3);
	t322 = t325 * t336 + t326 * t332;
	t321 = t322 * qJD(3);
	t333 = sin(qJ(2));
	t337 = cos(qJ(2));
	t348 = t322 * t333 + t340 * t337;
	t350 = t348 * qJD(2) + t320 * t337 + t333 * t321;
	t341 = t322 * t337 - t340 * t333;
	t329 = pkin(17) + pkin(18);
	t327 = sin(t329);
	t328 = cos(t329);
	t317 = -t327 * t341 - t348 * t328;
	t315 = 0.1e1 / t317 ^ 2;
	t316 = t327 * t348 - t341 * t328;
	t318 = -t341 * qJD(2) + t333 * t320 - t321 * t337;
	t310 = qJD(2) + qJD(3) + ((-t318 * t327 + t350 * t328) / t317 - (t318 * t328 + t327 * t350) * t316 * t315) / (t315 * t316 ^ 2 + 0.1e1);
	t313 = qJ(2) + qJ(3) + atan2(t316, t317);
	t346 = t310 * cos(t313);
	t311 = sin(t313);
	t342 = qJD(1) * t311;
	t338 = cos(qJ(1));
	t334 = sin(qJ(1));
	t1 = [0, 0, 0, t334 * t342 - t338 * t346; 0, 0, 0, -t334 * t346 - t338 * t342; 0, 0, 0, -t310 * t311;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobigD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:54
	% EndTime: 2020-05-07 02:13:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t128 = qJD(1) * cos(qJ(1));
	t127 = qJD(1) * sin(qJ(1));
	t1 = [0, t128, t128, 0; 0, t127, t127, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
end