% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh2m1OL
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh2m1OL_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1OL_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_jacobigD_rot_sym_varpar: pkin has to be [6x1] (double)');
JgD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, -qJD(1) * cos(qJ(1)), 0, 0, 0; 0, -qJD(1) * sin(qJ(1)), 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (4->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [0, -t30, -t30, 0, 0; 0, -t31, -t31, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (6->2), div. (0->0), fcn. (6->2), ass. (0->3)
	t38 = qJD(1) * sin(qJ(1));
	t37 = qJD(1) * cos(qJ(1));
	t1 = [0, -t37, -t37, -t37, 0; 0, -t38, -t38, -t38, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (25->11), mult. (15->8), div. (0->0), fcn. (15->4), ass. (0->9)
	t118 = qJD(2) + qJD(3) + qJD(4);
	t119 = qJ(2) + qJ(3) + qJ(4);
	t124 = cos(t119) * t118;
	t120 = sin(qJ(1));
	t123 = qJD(1) * t120;
	t121 = cos(qJ(1));
	t122 = qJD(1) * t121;
	t116 = sin(t119);
	t1 = [0, -t122, -t122, -t122, -t116 * t123 + t121 * t124; 0, -t123, -t123, -t123, t116 * t122 + t120 * t124; 0, 0, 0, 0, -t118 * t116;];
	JgD_rot = t1;
end