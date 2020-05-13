% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh1m2DE2
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
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JgD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh1m2DE2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_jacobigD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_jacobigD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_jacobigD_rot_sym_varpar: pkin has to be [22x1] (double)');
JgD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t29 = qJD(1) * cos(qJ(1));
	t28 = qJD(1) * sin(qJ(1));
	t1 = [0, t29, t29, 0; 0, t28, t28, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:31
	% EndTime: 2020-05-02 21:08:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:34
	% EndTime: 2020-05-02 21:08:34
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (752->26), mult. (1377->43), div. (18->3), fcn. (1980->15), ass. (0->33)
	t304 = pkin(22) + pkin(21);
	t302 = sin(t304);
	t303 = cos(t304);
	t305 = cos(pkin(20));
	t309 = sin(pkin(18));
	t324 = sin(pkin(20));
	t325 = cos(pkin(18));
	t299 = t309 * t305 - t325 * t324;
	t300 = t325 * t305 + t309 * t324;
	t306 = sin(qJ(3));
	t310 = cos(qJ(3));
	t295 = t299 * t310 - t306 * t300;
	t297 = t306 * t299 + t300 * t310;
	t307 = sin(qJ(2));
	t311 = cos(qJ(2));
	t327 = t307 * t295 + t297 * t311;
	t328 = -t295 * t311 + t307 * t297;
	t329 = -t328 * t302 + t327 * t303;
	t291 = t302 * t327 + t328 * t303;
	t290 = 0.1e1 / t329 ^ 2;
	t293 = t297 * qJD(3);
	t294 = t295 * qJD(3);
	t315 = qJD(3) * t311;
	t319 = t307 * t294;
	t320 = t307 * t293;
	t285 = qJD(2) + qJD(3) + (-((t297 * t315 + t319) * t303 + t302 * (t294 * t311 - t320) + t329 * qJD(2)) / t329 - ((-t295 * t315 + t320) * t303 - (-t293 * t311 - t319) * t302 + t291 * qJD(2)) * t291 * t290) / (t291 ^ 2 * t290 + 0.1e1);
	t288 = qJ(2) + qJ(3) + atan2(t291, -t329);
	t323 = t285 * cos(t288);
	t286 = sin(t288);
	t316 = qJD(1) * t286;
	t312 = cos(qJ(1));
	t308 = sin(qJ(1));
	t1 = [0, 0, 0, -t308 * t316 + t312 * t323; 0, 0, 0, t308 * t323 + t312 * t316; 0, 0, 0, t285 * t286;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobigD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t45 = qJD(1) * cos(qJ(1));
	t44 = qJD(1) * sin(qJ(1));
	t1 = [0, t45, t45, 0; 0, t44, t44, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobigD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobigD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t138 = qJD(1) * cos(qJ(1));
	t137 = qJD(1) * sin(qJ(1));
	t1 = [0, t138, t138, 0; 0, t137, t137, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
end