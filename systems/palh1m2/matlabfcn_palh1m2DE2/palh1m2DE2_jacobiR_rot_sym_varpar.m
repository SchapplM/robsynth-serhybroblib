% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh1m2DE2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_jacobiR_rot_sym_varpar: pkin has to be [22x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:27
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0; -t9, 0, 0, 0; -t8, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t12 = cos(qJ(1));
	t9 = sin(qJ(2));
	t15 = t12 * t9;
	t10 = sin(qJ(1));
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t8 = t10 * t9;
	t1 = [t8, -t13, 0, 0; -t15, -t14, 0, 0; 0, -t9, 0, 0; t14, t15, 0, 0; -t13, t8, 0, 0; 0, -t11, 0, 0; t12, 0, 0, 0; t10, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.02s
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
	t1 = [-t26, -t25, -t25, 0; t24, -t27, -t27, 0; 0, t20, t20, 0; t27, -t24, -t24, 0; -t25, -t26, -t26, 0; 0, -t19, -t19, 0; t23, 0, 0, 0; t22, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:31
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (2472->13), mult. (4296->20), div. (72->0), fcn. (6486->15), ass. (0->23)
	t101 = sin(pkin(20));
	t110 = cos(pkin(18));
	t93 = cos(pkin(20));
	t97 = sin(pkin(18));
	t87 = -t110 * t101 + t97 * t93;
	t88 = t97 * t101 + t110 * t93;
	t94 = sin(qJ(3));
	t98 = cos(qJ(3));
	t112 = t94 * t87 + t88 * t98;
	t82 = t87 * t98 - t94 * t88;
	t95 = sin(qJ(2));
	t99 = cos(qJ(2));
	t111 = -t95 * t112 + t82 * t99;
	t113 = t112 * t99 + t95 * t82;
	t100 = cos(qJ(1));
	t96 = sin(qJ(1));
	t92 = pkin(22) + pkin(21);
	t91 = cos(t92);
	t90 = sin(t92);
	t70 = qJ(2) + qJ(3) + atan2(-t111 * t91 + t113 * t90, -t111 * t90 - t113 * t91);
	t69 = cos(t70);
	t68 = sin(t70);
	t1 = [-t96 * t69, 0, 0, 0; t100 * t69, 0, 0, 0; 0, 0, 0, 0; t96 * t68, 0, 0, 0; -t100 * t68, 0, 0, 0; 0, 0, 0, 0; t100, 0, 0, 0; t96, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:33
	% EndTime: 2020-05-02 21:08:34
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (3923->20), mult. (6838->28), div. (108->0), fcn. (10332->17), ass. (0->33)
	t246 = cos(pkin(20));
	t251 = sin(pkin(18));
	t272 = sin(pkin(20));
	t273 = cos(pkin(18));
	t240 = t246 * t251 - t272 * t273;
	t241 = t246 * t273 + t251 * t272;
	t248 = sin(qJ(3));
	t253 = cos(qJ(3));
	t235 = t240 * t253 - t241 * t248;
	t249 = sin(qJ(2));
	t254 = cos(qJ(2));
	t275 = t240 * t248 + t241 * t253;
	t280 = t235 * t249 + t254 * t275;
	t274 = -t235 * t254 + t249 * t275;
	t247 = sin(qJ(4));
	t250 = sin(qJ(1));
	t263 = t250 * t247;
	t252 = cos(qJ(4));
	t262 = t250 * t252;
	t255 = cos(qJ(1));
	t261 = t255 * t247;
	t260 = t255 * t252;
	t245 = pkin(22) + pkin(21);
	t244 = cos(t245);
	t243 = sin(t245);
	t223 = qJ(2) + qJ(3) + atan2(t280 * t243 + t274 * t244, t274 * t243 - t244 * t280);
	t222 = cos(t223);
	t221 = sin(t223);
	t220 = t222 * t260 + t263;
	t219 = -t222 * t261 + t262;
	t218 = -t222 * t262 + t261;
	t217 = t222 * t263 + t260;
	t1 = [t218, 0, 0, t219; t220, 0, 0, -t217; 0, 0, 0, -t221 * t247; t217, 0, 0, -t220; t219, 0, 0, t218; 0, 0, 0, -t221 * t252; -t250 * t221, 0, 0, 0; t255 * t221, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->8), mult. (68->12), div. (0->0), fcn. (110->8), ass. (0->17)
	t28 = sin(pkin(18));
	t29 = sin(pkin(17));
	t32 = cos(pkin(18));
	t33 = cos(pkin(17));
	t24 = t28 * t33 - t32 * t29;
	t25 = t29 * t28 + t33 * t32;
	t26 = sin(qJ(2));
	t30 = cos(qJ(2));
	t20 = -t26 * t24 - t25 * t30;
	t27 = sin(qJ(1));
	t37 = t20 * t27;
	t21 = t24 * t30 - t26 * t25;
	t36 = t21 * t27;
	t31 = cos(qJ(1));
	t35 = t31 * t20;
	t34 = t31 * t21;
	t1 = [-t36, t35, 0, 0; t34, t37, 0, 0; 0, t21, 0, 0; -t37, -t34, 0, 0; t35, -t36, 0, 0; 0, t20, 0, 0; t31, 0, 0, 0; t27, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->7), mult. (8->6), div. (0->0), fcn. (10->6), ass. (0->8)
	t21 = pkin(18) - pkin(22);
	t20 = -qJ(1) + t21;
	t19 = qJ(1) + t21;
	t18 = cos(t19);
	t17 = sin(t20);
	t16 = -cos(t20) / 0.2e1;
	t15 = sin(t19) / 0.2e1;
	t1 = [-t17 / 0.2e1 + t15, 0, 0, 0; t16 - t18 / 0.2e1, 0, 0, 0; 0, 0, 0, 0; t18 / 0.2e1 + t16, 0, 0, 0; t15 + t17 / 0.2e1, 0, 0, 0; 0, 0, 0, 0; cos(qJ(1)), 0, 0, 0; sin(qJ(1)), 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (118->13), mult. (202->8), div. (30->0), fcn. (342->9), ass. (0->14)
	t43 = sin(qJ(3));
	t35 = cos(pkin(19));
	t37 = cos(qJ(3));
	t39 = sin(pkin(19));
	t28 = qJ(2) + atan2(-t37 * t35 + t43 * t39, t43 * t35 + t37 * t39);
	t26 = sin(t28);
	t36 = sin(qJ(1));
	t25 = t36 * t26;
	t27 = cos(t28);
	t42 = t36 * t27;
	t38 = cos(qJ(1));
	t41 = t38 * t26;
	t40 = t38 * t27;
	t1 = [t25, -t40, -t40, 0; -t41, -t42, -t42, 0; 0, -t26, -t26, 0; t42, t41, t41, 0; -t40, t25, t25, 0; 0, -t27, -t27, 0; t38, 0, 0, 0; t36, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiR_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (234->11), mult. (398->8), div. (66->0), fcn. (678->10), ass. (0->16)
	t67 = sin(qJ(3));
	t57 = cos(pkin(19));
	t59 = cos(qJ(3));
	t62 = sin(pkin(19));
	t52 = t67 * t57 + t59 * t62;
	t54 = t59 * t57 - t67 * t62;
	t46 = qJ(2) + atan2(-t54, t52) + atan2(-t54, -t52);
	t44 = sin(t46);
	t58 = sin(qJ(1));
	t65 = t58 * t44;
	t45 = cos(t46);
	t64 = t58 * t45;
	t60 = cos(qJ(1));
	t63 = t60 * t44;
	t43 = t60 * t45;
	t1 = [-t65, t43, 0, 0; t63, t64, 0, 0; 0, t44, 0, 0; -t64, -t63, 0, 0; t43, -t65, 0, 0; 0, t45, 0, 0; t60, 0, 0, 0; t58, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiR_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (1644->18), mult. (2608->24), div. (108->0), fcn. (4002->16), ass. (0->27)
	t108 = sin(qJ(2));
	t110 = cos(qJ(2));
	t122 = sin(qJ(3));
	t123 = sin(pkin(18));
	t124 = cos(qJ(3));
	t125 = cos(pkin(18));
	t98 = -t122 * t123 - t124 * t125;
	t99 = -t122 * t125 + t124 * t123;
	t113 = t108 * t98 + t110 * t99;
	t109 = sin(qJ(1));
	t106 = pkin(22) + pkin(21);
	t104 = sin(t106);
	t105 = cos(t106);
	t112 = t108 * t99 - t98 * t110;
	t107 = sin(pkin(22));
	t114 = cos(pkin(22));
	t96 = t123 * t107 + t114 * t125;
	t97 = t125 * t107 - t123 * t114;
	t77 = -qJ(2) - atan2(-t97 * t108 + t96 * t110, t108 * t96 + t97 * t110) + pkin(21) - atan2(t104 * t112 - t105 * t113, t104 * t113 + t112 * t105);
	t75 = sin(t77);
	t118 = t109 * t75;
	t76 = cos(t77);
	t117 = t109 * t76;
	t111 = cos(qJ(1));
	t116 = t111 * t75;
	t115 = t111 * t76;
	t1 = [t118, t115, t115, 0; -t116, t117, t117, 0; 0, -t75, -t75, 0; -t117, t116, t116, 0; t115, t118, t118, 0; 0, t76, t76, 0; t111, 0, 0, 0; t109, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end