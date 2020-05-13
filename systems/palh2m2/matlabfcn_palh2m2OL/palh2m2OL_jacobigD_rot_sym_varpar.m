% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh2m2OL
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh2m2OL_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2OL_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_jacobigD_rot_sym_varpar: pkin has to be [5x1] (double)');
JgD_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t28 = qJD(1) * cos(qJ(1));
	t27 = qJD(1) * sin(qJ(1));
	t1 = [0, t28, t28, 0, 0, 0; 0, t27, t27, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (6->2), div. (0->0), fcn. (6->2), ass. (0->3)
	t34 = qJD(1) * cos(qJ(1));
	t33 = qJD(1) * sin(qJ(1));
	t1 = [0, t34, t34, t34, 0, 0; 0, t33, t33, t33, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->3)
	t40 = qJD(1) * cos(qJ(1));
	t39 = qJD(1) * sin(qJ(1));
	t1 = [0, t40, t40, t40, t40, 0; 0, t39, t39, t39, t39, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:58
	% EndTime: 2020-05-03 06:34:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->5), mult. (17->8), div. (0->0), fcn. (17->4), ass. (0->9)
	t111 = qJD(2) + qJD(3) + qJD(4) + qJD(5);
	t112 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
	t117 = cos(t112) * t111;
	t115 = sin(qJ(1));
	t113 = qJD(1) * t115;
	t116 = cos(qJ(1));
	t114 = qJD(1) * t116;
	t109 = sin(t112);
	t1 = [0, t114, t114, t114, t114, t109 * t113 - t116 * t117; 0, t113, t113, t113, t113, -t109 * t114 - t115 * t117; 0, 0, 0, 0, 0, -t111 * t109;];
	JgD_rot = t1;
end