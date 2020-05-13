% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh1m2OL
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% Jg_rot [3x13]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh1m2OL_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),uint8(0),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_jacobig_rot_sym_varpar: qJ has to be [13x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2OL_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_jacobig_rot_sym_varpar: pkin has to be [20x1] (double)');
Jg_rot=NaN(3,13);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:40
	% EndTime: 2020-05-02 23:30:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:40
	% EndTime: 2020-05-02 23:30:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:40
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -cos(qJ(1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t18 = cos(qJ(1));
	t17 = sin(qJ(1));
	t1 = [0, t17, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t18, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (6->2), ass. (0->3)
	t24 = cos(qJ(1));
	t23 = sin(qJ(1));
	t1 = [0, t23, t23, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t24, -t24, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->5), mult. (2->2), div. (0->0), fcn. (11->4), ass. (0->5)
	t94 = cos(qJ(1));
	t93 = sin(qJ(1));
	t92 = qJ(2) + qJ(3) + qJ(4);
	t91 = sin(t92);
	t1 = [0, t93, t93, t93, t94 * t91, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t94, -t94, -t94, t93 * t91, 0, 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, -cos(t92), 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, sin(qJ(1)), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -cos(qJ(1)), 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:40
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t19 = cos(qJ(1));
	t18 = sin(qJ(1));
	t1 = [0, t18, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0; 0, -t19, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:40
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [0, t11, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0; 0, -t12, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobig_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (6->2), ass. (0->3)
	t24 = cos(qJ(1));
	t23 = sin(qJ(1));
	t1 = [0, t23, 0, 0, 0, 0, 0, t23, t23, 0, 0, 0, 0; 0, -t24, 0, 0, 0, 0, 0, -t24, -t24, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobig_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (6->2), ass. (0->3)
	t30 = cos(qJ(1));
	t29 = sin(qJ(1));
	t1 = [0, t29, 0, 0, 0, 0, t29, 0, 0, t29, 0, 0, 0; 0, -t30, 0, 0, 0, 0, -t30, 0, 0, -t30, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
end