% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh1m2TE
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
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh1m2TE_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_jacobigD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_jacobigD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2TE_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_jacobigD_rot_sym_varpar: pkin has to be [22x1] (double)');
JgD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = qJD(1) * cos(qJ(1));
	t35 = qJD(1) * sin(qJ(1));
	t1 = [0, t36, t36, 0; 0, t35, t35, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->5), mult. (16->9), div. (0->0), fcn. (22->8), ass. (0->7)
	t106 = pkin(22) + pkin(21);
	t107 = sin(pkin(20));
	t108 = cos(pkin(20));
	t109 = sin(pkin(18));
	t110 = cos(pkin(18));
	t111 = qJD(1) * ((t110 * t107 - t109 * t108) * cos(t106) + sin(t106) * (t109 * t107 + t110 * t108));
	t1 = [0, 0, 0, -sin(qJ(1)) * t111; 0, 0, 0, cos(qJ(1)) * t111; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobigD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t75 = qJD(1) * cos(qJ(1));
	t74 = qJD(1) * sin(qJ(1));
	t1 = [0, t75, t75, 0; 0, t74, t74, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobigD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobigD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t34 = qJD(1) * cos(qJ(1));
	t33 = qJD(1) * sin(qJ(1));
	t1 = [0, t34, t34, 0; 0, t33, t33, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
end