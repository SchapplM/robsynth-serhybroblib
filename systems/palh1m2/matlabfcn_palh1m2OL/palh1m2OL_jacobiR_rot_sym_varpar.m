% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
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
% JR_rot [9x13]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh1m2OL_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),uint8(0),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_jacobiR_rot_sym_varpar: qJ has to be [13x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2OL_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_jacobiR_rot_sym_varpar: pkin has to be [20x1] (double)');
JR_rot=NaN(9,13);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:40
	% EndTime: 2020-05-02 23:30:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:40
	% EndTime: 2020-05-02 23:30:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t12 = cos(qJ(1));
	t9 = sin(qJ(2));
	t15 = t12 * t9;
	t10 = sin(qJ(1));
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t8 = t10 * t9;
	t1 = [t8, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t13, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t21 = qJ(2) + qJ(3);
	t19 = sin(t21);
	t22 = sin(qJ(1));
	t27 = t22 * t19;
	t20 = cos(t21);
	t26 = t22 * t20;
	t23 = cos(qJ(1));
	t25 = t23 * t19;
	t24 = t23 * t20;
	t1 = [-t26, -t25, -t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t24, -t27, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t20, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t27, -t24, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t25, -t26, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t19, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->18), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t27 = qJ(2) + qJ(3) + qJ(4);
	t25 = sin(t27);
	t28 = sin(qJ(1));
	t33 = t28 * t25;
	t26 = cos(t27);
	t32 = t28 * t26;
	t29 = cos(qJ(1));
	t31 = t29 * t25;
	t30 = t29 * t26;
	t1 = [-t32, -t31, -t31, -t31, 0, 0, 0, 0, 0, 0, 0, 0, 0; t30, -t33, -t33, -t33, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t26, t26, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0; t33, -t30, -t30, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t31, -t32, -t32, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t25, -t25, -t25, 0, 0, 0, 0, 0, 0, 0, 0, 0; t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:42
	% EndTime: 2020-05-02 23:30:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (98->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
	t106 = qJ(2) + qJ(3) + qJ(4);
	t105 = cos(t106);
	t107 = sin(qJ(5));
	t117 = t105 * t107;
	t108 = sin(qJ(1));
	t116 = t108 * t107;
	t109 = cos(qJ(5));
	t115 = t108 * t109;
	t110 = cos(qJ(1));
	t114 = t110 * t107;
	t113 = t110 * t109;
	t104 = sin(t106);
	t112 = t104 * t115;
	t111 = t104 * t113;
	t103 = t110 * t105;
	t102 = t105 * t109;
	t101 = t108 * t105;
	t100 = t104 * t114;
	t99 = t104 * t116;
	t98 = t105 * t113 + t116;
	t97 = -t105 * t114 + t115;
	t96 = -t105 * t115 + t114;
	t95 = t105 * t116 + t113;
	t1 = [t96, -t111, -t111, -t111, t97, 0, 0, 0, 0, 0, 0, 0, 0; t98, -t112, -t112, -t112, -t95, 0, 0, 0, 0, 0, 0, 0, 0; 0, t102, t102, t102, -t104 * t107, 0, 0, 0, 0, 0, 0, 0, 0; t95, t100, t100, t100, -t98, 0, 0, 0, 0, 0, 0, 0, 0; t97, t99, t99, t99, t96, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t117, -t117, -t117, -t104 * t109, 0, 0, 0, 0, 0, 0, 0, 0; -t108 * t104, t103, t103, t103, 0, 0, 0, 0, 0, 0, 0, 0, 0; t110 * t104, t101, t101, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t104, t104, t104, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(6));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(6));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0; t12, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t24 = qJ(2) + qJ(7);
	t23 = cos(t24);
	t25 = sin(qJ(1));
	t28 = t25 * t23;
	t26 = cos(qJ(1));
	t27 = t26 * t23;
	t22 = sin(t24);
	t21 = t26 * t22;
	t20 = t25 * t22;
	t1 = [t20, -t27, 0, 0, 0, 0, -t27, 0, 0, 0, 0, 0, 0; -t21, -t28, 0, 0, 0, 0, -t28, 0, 0, 0, 0, 0, 0; 0, -t22, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0; t28, t21, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0; -t27, t20, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0; 0, -t23, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0; t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t17 = qJ(2) + qJ(8);
	t16 = cos(t17);
	t18 = sin(qJ(1));
	t21 = t18 * t16;
	t19 = cos(qJ(1));
	t20 = t19 * t16;
	t15 = sin(t17);
	t14 = t19 * t15;
	t13 = t18 * t15;
	t1 = [t13, -t20, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0; -t14, -t21, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0; 0, -t15, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0; t21, t14, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0; -t20, t13, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0; 0, -t16, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiR_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (52->9), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t29 = qJ(2) + qJ(8) + qJ(9);
	t27 = sin(t29);
	t30 = sin(qJ(1));
	t33 = t30 * t27;
	t31 = cos(qJ(1));
	t32 = t31 * t27;
	t28 = cos(t29);
	t26 = t31 * t28;
	t25 = t30 * t28;
	t1 = [-t33, t26, 0, 0, 0, 0, 0, t26, t26, 0, 0, 0, 0; t32, t25, 0, 0, 0, 0, 0, t25, t25, 0, 0, 0, 0; 0, t27, 0, 0, 0, 0, 0, t27, t27, 0, 0, 0, 0; -t25, -t32, 0, 0, 0, 0, 0, -t32, -t32, 0, 0, 0, 0; t26, -t33, 0, 0, 0, 0, 0, -t33, -t33, 0, 0, 0, 0; 0, t28, 0, 0, 0, 0, 0, t28, t28, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiR_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:30:41
	% EndTime: 2020-05-02 23:30:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (71->6), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t37 = -qJ(2) - qJ(7) + pkin(19) - qJ(10);
	t36 = cos(t37);
	t35 = sin(t37);
	t34 = t39 * t36;
	t33 = t39 * t35;
	t32 = t38 * t36;
	t31 = t38 * t35;
	t1 = [t31, t34, 0, 0, 0, 0, t34, 0, 0, t34, 0, 0, 0; -t33, t32, 0, 0, 0, 0, t32, 0, 0, t32, 0, 0, 0; 0, -t35, 0, 0, 0, 0, -t35, 0, 0, -t35, 0, 0, 0; -t32, t33, 0, 0, 0, 0, t33, 0, 0, t33, 0, 0, 0; t34, t31, 0, 0, 0, 0, t31, 0, 0, t31, 0, 0, 0; 0, t36, 0, 0, 0, 0, t36, 0, 0, t36, 0, 0, 0; t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
end