% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% mg10hlDE1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in mg10hlDE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 12:53
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = mg10hlDE1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlDE1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE1_jacobia_rot_sym_varpar: pkin has to be [17x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:05
	% EndTime: 2020-04-11 12:36:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Ja_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:05
	% EndTime: 2020-04-11 12:36:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (6->3), div. (5->2), fcn. (6->2), ass. (0->24)
	unknown=NaN(3,6);
	t1 = sin(qJ(1));
	t2 = t1 ^ 2;
	t3 = cos(qJ(1));
	t4 = t3 ^ 2;
	t6 = t2 / t4;
	t8 = 0.1e1 / (0.1e1 + t6);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(3,1) = (t6 * t8 + t8);
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	Ja_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:05
	% EndTime: 2020-04-11 12:36:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:08
	% EndTime: 2020-04-11 12:36:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:10
	% EndTime: 2020-04-11 12:36:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:16
	% EndTime: 2020-04-11 12:36:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:22
	% EndTime: 2020-04-11 12:36:49
	% DurationCPUTime: 10.50s
	% Computational Cost: add. (563083->201), mult. (746406->475), div. (56329->32), fcn. (420980->26), ass. (0->227)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t2 * t1;
	t4 = cos(pkin(15));
	t5 = cos(qJ(2));
	t7 = sin(pkin(15));
	t9 = -t7 * t2 + t4 * t5;
	t10 = t9 * pkin(2);
	t13 = t4 * t2 + t7 * t5;
	t16 = 0.2e1 * pkin(2) * pkin(1) * t13;
	t17 = pkin(1) ^ 2;
	t21 = -t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t25 = -t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t27 = sqrt(-t25 * t21);
	t28 = t27 * t10;
	t30 = -t13 * pkin(2) + pkin(1);
	t31 = pkin(2) ^ 2;
	t32 = pkin(5) ^ 2;
	t33 = pkin(4) ^ 2;
	t34 = -t16 + t17 + t31 - t32 + t33;
	t36 = t34 * t30 - t28;
	t37 = t36 * t4;
	t38 = 0.1e1 / pkin(4);
	t39 = -t16 + t17 + t31;
	t40 = 0.1e1 / t39;
	t43 = t34 * t10;
	t44 = t27 * t30 + t43;
	t45 = t44 ^ 2;
	t46 = 0.1e1 / t33;
	t47 = t46 * t45;
	t48 = t39 ^ 2;
	t49 = 0.1e1 / t48;
	t51 = t36 ^ 2;
	t52 = t46 * t51;
	t54 = t49 * t47 + t49 * t52;
	t55 = sqrt(t54);
	t56 = 0.1e1 / t55;
	t57 = t56 * t40 * t38;
	t59 = t44 * t7;
	t61 = -t57 * t37 + t57 * t59;
	t62 = t61 * t3;
	t63 = t5 * t1;
	t64 = t36 * t7;
	t66 = t44 * t4;
	t68 = -t57 * t64 - t57 * t66;
	t70 = -t68 * t63 - t62;
	t71 = cos(pkin(17));
	t72 = 1 / pkin(7);
	t73 = qJ(6) + pkin(8);
	t75 = 0.1e1 / t73 * t72;
	t76 = (pkin(6) - pkin(7) - pkin(8) - qJ(6));
	t77 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t78 = t77 * t76;
	t79 = (pkin(6) + pkin(7) - pkin(8) - qJ(6));
	t80 = (pkin(6) + pkin(7) + pkin(8) + qJ(6));
	t81 = t80 * t79;
	t83 = sqrt(-t81 * t78);
	t85 = (pkin(6) ^ 2);
	t86 = (pkin(7) ^ 2);
	t87 = (pkin(8) ^ 2);
	t89 = 2 * pkin(8) * qJ(6);
	t90 = qJ(6) ^ 2;
	t91 = t85 - t86 - t87 - t89 - t90;
	t93 = atan2(t83 * t75, t91 * t75);
	t94 = t93 + pkin(16);
	t95 = cos(t94);
	t97 = sin(pkin(17));
	t98 = sin(t94);
	t100 = -t95 * t71 - t98 * t97;
	t103 = t61 * t63;
	t104 = t68 * t3 - t103;
	t107 = t98 * t71 - t95 * t97;
	t109 = t100 * t70 + t107 * t104;
	t110 = 0.1e1 / pkin(6);
	t111 = t73 ^ 2;
	t112 = 1 / t111;
	t113 = (t112 * t110);
	t114 = t85 - t86 + t87 + t89 + t90;
	t115 = t72 * t114;
	t116 = t91 * t115;
	t118 = (t76 * t113);
	t119 = (t79 * t77);
	t121 = t72 * t80 * t119;
	t123 = t116 * t113 - t121 * t118;
	t124 = cos(pkin(16));
	t126 = t72 * t83;
	t127 = t91 * t126;
	t129 = t83 * t115;
	t131 = -t127 * t113 + t129 * t113;
	t132 = sin(pkin(16));
	t134 = -t124 * t123 / 0.4e1 + t132 * t131 / 0.4e1;
	t135 = t134 * t109;
	t138 = t124 * t131 / 0.4e1 + t132 * t123 / 0.4e1;
	t139 = t138 ^ 2;
	t140 = t134 ^ 2;
	t141 = t139 + t140;
	t142 = sqrt(t141);
	t143 = 0.1e1 / t142;
	t147 = t100 * t104 - t107 * t70;
	t148 = t138 * t147;
	t150 = t143 * t135 + t143 * t148;
	t151 = t61 * t5;
	t153 = -t68 * t2 + t151;
	t156 = t61 * t2;
	t157 = -t68 * t5 - t156;
	t159 = t100 * t153 + t107 * t157;
	t160 = t134 * t159;
	t164 = t100 * t157 - t107 * t153;
	t165 = t138 * t164;
	t167 = -t143 * t160 - t143 * t165;
	t168 = 0.1e1 / t167;
	t169 = t168 * t150;
	t170 = sin(qJ(1));
	t171 = t2 * t170;
	t172 = t61 * t171;
	t173 = t5 * t170;
	t175 = -t68 * t173 - t172;
	t178 = t61 * t173;
	t179 = t68 * t171 - t178;
	t181 = t100 * t175 + t107 * t179;
	t182 = t134 * t181;
	t186 = t100 * t179 - t107 * t175;
	t187 = t138 * t186;
	t189 = t143 * t182 + t143 * t187;
	t190 = t189 ^ 2;
	t191 = t167 ^ 2;
	t192 = 0.1e1 / t191;
	t195 = 0.1e1 / (t192 * t190 + 0.1e1);
	t197 = -t13 * pkin(2);
	t199 = 0.1e1 / t27;
	t200 = pkin(1) * t9;
	t204 = pkin(1) * pkin(2);
	t206 = t25 * pkin(2) * t200 + t204 * t9 * t21;
	t213 = -t206 * t199 * t10 - 0.2e1 * t204 * t9 * t30 - t27 * t197 - t43;
	t216 = t49 * t38;
	t219 = t204 * t9 * t56;
	t231 = t9 ^ 2;
	t235 = -0.2e1 * pkin(1) * t231 * t31 + t206 * t199 * t30 + t34 * t197 - t28;
	t240 = 0.1e1 / t48 / t39;
	t242 = pkin(2) * t200;
	t253 = (0.2e1 * t213 * t49 * t46 * t36 + 0.2e1 * t235 * t49 * t46 * t44 + 0.4e1 * t242 * t240 * t47 + 0.4e1 * t242 * t240 * t52) / t55 / t54 * t40;
	t264 = -t57 * t213 * t4 - 0.2e1 * t219 * t216 * t37 + t253 * t38 * t37 / 0.2e1 + t57 * t235 * t7 + 0.2e1 * t219 * t216 * t59 - t253 * t38 * t59 / 0.2e1;
	t283 = -t57 * t213 * t7 - 0.2e1 * t219 * t216 * t64 + t253 * t38 * t64 / 0.2e1 - t57 * t235 * t4 - 0.2e1 * t219 * t216 * t66 + t253 * t38 * t66 / 0.2e1;
	t285 = -t264 * t171 + t68 * t171 - t283 * t173 - t178;
	t290 = t283 * t171 - t264 * t173 + t68 * t173 + t172;
	t300 = t143 * t134 * (t100 * t285 + t107 * t290) + t143 * t138 * (t100 * t290 - t107 * t285);
	t306 = -t283 * t2 + t264 * t5 - t68 * t5 - t156;
	t311 = -t264 * t2 + t68 * t2 - t283 * t5 - t151;
	t321 = -t143 * t134 * (t100 * t306 + t107 * t311) - t143 * t138 * (t100 * t311 - t107 * t306);
	t323 = t195 * t192;
	t325 = t195 * t168 * t300 - t323 * t189 * t321;
	t326 = t112 * t72;
	t328 = 0.1e1 / t83;
	t334 = -(t80 * t79 * t76) + (t80 * t119) - t79 * t78 + t80 * t78;
	t342 = t91 ^ 2;
	t343 = 1 / t342;
	t347 = 0.1e1 / (-t343 * t81 * t78 + 0.1e1);
	t358 = t347 / t91 * t73 * pkin(7) * (-t83 * t326 + t334 * t328 * t75 / 0.2e1) - t347 * t343 * t83 * t73 * pkin(7) * (-(t91 * t326) - 0.2e1 * t73 * t75);
	t359 = t358 * t71;
	t361 = t358 * t97;
	t363 = t98 * t359 - t95 * t361;
	t367 = t95 * t359 + t98 * t361;
	t374 = 0.1e1 / t111 / t73 * t110;
	t377 = 0.2e1 * t72 * t73;
	t388 = t72 * t81;
	t400 = -t116 * t374 / 0.2e1 + t91 * t377 * t113 / 0.4e1 - t73 * t115 * t113 / 0.2e1 + t121 * t76 * t374 / 0.2e1 + t388 * t77 * t113 / 0.4e1 - t388 * t118 / 0.4e1 + t72 * t80 * t77 * t118 / 0.4e1 - (t72 * t119 * t118) / 0.4e1;
	t422 = t127 * t374 / 0.2e1 - t334 * t91 * t72 * t328 * t113 / 0.8e1 + t73 * t126 * t113 / 0.2e1 - t129 * t374 / 0.2e1 + t83 * t377 * t113 / 0.4e1 + t334 * t328 * t72 * t114 * t113 / 0.8e1;
	t424 = -t124 * t400 + t132 * t422;
	t431 = t124 * t422 + t132 * t400;
	t435 = 0.2e1 * (t424 * t134 + t431 * t138) / t142 / t141;
	t447 = t143 * t134 * (t363 * t175 + t367 * t179) + t143 * t424 * t181 - t435 * t182 / 0.2e1 + t143 * t138 * (-t367 * t175 + t363 * t179) + t143 * t431 * t186 - t435 * t187 / 0.2e1;
	t468 = -t143 * t134 * (t363 * t153 + t367 * t157) - t143 * t424 * t159 + t435 * t160 / 0.2e1 - t143 * t138 * (-t367 * t153 + t363 * t157) - t143 * t431 * t164 + t435 * t165 / 0.2e1;
	t471 = t195 * t168 * t447 - t323 * t189 * t468;
	t474 = -t100 * t175 - t107 * t179;
	t479 = -t100 * t179 + t107 * t175;
	t483 = atan2(t189, t167);
	t484 = cos(t483);
	t486 = sin(t483);
	t488 = t167 * t484 + t189 * t486;
	t489 = 0.1e1 / t488;
	t491 = t150 ^ 2;
	t492 = t488 ^ 2;
	t493 = 0.1e1 / t492;
	t496 = 0.1e1 / (t493 * t491 + 0.1e1);
	t506 = t496 * t493;
	t512 = -t264 * t3 - t283 * t63 + t68 * t3 - t103;
	t517 = -t264 * t63 + t283 * t3 + t68 * t63 + t62;
	t519 = t100 * t512 + t107 * t517;
	t524 = t100 * t517 - t107 * t512;
	t542 = t367 * t104 + t363 * t70;
	t551 = t363 * t104 - t367 * t70;
	t575 = -t143 * t134 * t479 + t143 * t138 * t474;
	t576 = sin(qJ(3));
	t578 = cos(qJ(3));
	t581 = t138 * t109;
	t583 = t134 * t147;
	t585 = t143 * t581 - t143 * t583;
	t588 = t576 * t170 + t578 * t585;
	t589 = 0.1e1 / t588;
	t593 = -t578 * t170 + t576 * t585;
	t594 = t593 ^ 2;
	t595 = t588 ^ 2;
	t596 = 0.1e1 / t595;
	t599 = 0.1e1 / (t596 * t594 + 0.1e1);
	t605 = t599 * t596;
	t612 = -t143 * t134 * t524 + t143 * t138 * t519;
	t614 = t599 * t589;
	t618 = t599 * t596 * t593;
	t636 = t143 * t138 * t542 + t143 * t431 * t109 - t435 * t581 / 0.2e1 - t143 * t134 * t551 - t143 * t424 * t147 + t435 * t583 / 0.2e1;
	unknown(1,1) = t195 * t169;
	unknown(1,2) = t325;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t471;
	unknown(2,1) = t496 * t489 * (-t143 * t134 * t474 - t143 * t138 * t479) + t506 * t150 * (t189 * t484 * t195 * t169 - t486 * t195 * t150 + t150 * t486);
	unknown(2,2) = t496 * t489 * (-t143 * t134 * t519 - t143 * t138 * t524) + t506 * t150 * (-t167 * t486 * t325 + t189 * t484 * t325 + t300 * t486 + t321 * t484);
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = t496 * t489 * (-t143 * t134 * t542 - t143 * t424 * t109 + t435 * t135 / 0.2e1 - t143 * t138 * t551 - t143 * t431 * t147 + t435 * t148 / 0.2e1) + t506 * t150 * (-t167 * t486 * t471 + t189 * t484 * t471 + t447 * t486 + t468 * t484);
	unknown(3,1) = t599 * t589 * (-t578 * t1 + t576 * t575) - t605 * t593 * (t576 * t1 + t578 * t575);
	unknown(3,2) = t614 * t576 * t612 - t618 * t578 * t612;
	unknown(3,3) = t605 * t593 ^ 2 + t599;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = t614 * t576 * t636 - t618 * t578 * t636;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:34
	% EndTime: 2020-04-11 12:37:05
	% DurationCPUTime: 14.68s
	% Computational Cost: add. (804940->219), mult. (1067488->527), div. (81602->34), fcn. (602271->28), ass. (0->247)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t2 * t1;
	t4 = cos(pkin(15));
	t5 = cos(qJ(2));
	t7 = sin(pkin(15));
	t9 = -t7 * t2 + t4 * t5;
	t10 = t9 * pkin(2);
	t13 = t4 * t2 + t7 * t5;
	t16 = 0.2e1 * pkin(2) * pkin(1) * t13;
	t17 = pkin(1) ^ 2;
	t21 = -t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t25 = -t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t27 = sqrt(-t25 * t21);
	t28 = t27 * t10;
	t30 = -t13 * pkin(2) + pkin(1);
	t31 = pkin(2) ^ 2;
	t32 = pkin(5) ^ 2;
	t33 = pkin(4) ^ 2;
	t34 = -t16 + t17 + t31 - t32 + t33;
	t36 = t34 * t30 - t28;
	t37 = t36 * t4;
	t38 = 0.1e1 / pkin(4);
	t39 = -t16 + t17 + t31;
	t40 = 0.1e1 / t39;
	t43 = t34 * t10;
	t44 = t27 * t30 + t43;
	t45 = t44 ^ 2;
	t46 = 0.1e1 / t33;
	t47 = t46 * t45;
	t48 = t39 ^ 2;
	t49 = 0.1e1 / t48;
	t51 = t36 ^ 2;
	t52 = t46 * t51;
	t54 = t49 * t47 + t49 * t52;
	t55 = sqrt(t54);
	t56 = 0.1e1 / t55;
	t57 = t56 * t40 * t38;
	t59 = t44 * t7;
	t61 = -t57 * t37 + t57 * t59;
	t62 = t61 * t3;
	t63 = t5 * t1;
	t64 = t36 * t7;
	t66 = t44 * t4;
	t68 = -t57 * t64 - t57 * t66;
	t70 = -t68 * t63 - t62;
	t71 = cos(pkin(17));
	t72 = 0.1e1 / pkin(7);
	t73 = qJ(6) + pkin(8);
	t75 = 0.1e1 / t73 * t72;
	t76 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t77 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t78 = t77 * t76;
	t79 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t80 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t81 = t80 * t79;
	t83 = sqrt(-t81 * t78);
	t85 = (pkin(6) ^ 2);
	t86 = (pkin(7) ^ 2);
	t87 = (pkin(8) ^ 2);
	t89 = 2 * pkin(8) * qJ(6);
	t90 = qJ(6) ^ 2;
	t91 = t85 - t86 - t87 - t89 - t90;
	t93 = atan2(t83 * t75, t91 * t75);
	t94 = t93 + pkin(16);
	t95 = cos(t94);
	t97 = sin(pkin(17));
	t98 = sin(t94);
	t100 = -t95 * t71 - t98 * t97;
	t103 = t61 * t63;
	t104 = t68 * t3 - t103;
	t107 = t98 * t71 - t95 * t97;
	t109 = t100 * t70 + t107 * t104;
	t110 = 0.1e1 / pkin(6);
	t111 = t73 ^ 2;
	t112 = 0.1e1 / t111;
	t113 = t112 * t110;
	t114 = t72 * t83;
	t115 = t91 * t114;
	t117 = t85 - t86 + t87 + t89 + t90;
	t118 = t72 * t117;
	t119 = t83 * t118;
	t121 = -t115 * t113 + t119 * t113;
	t122 = cos(pkin(16));
	t124 = t91 * t118;
	t126 = t76 * t113;
	t127 = t79 * t77;
	t129 = t72 * t80 * t127;
	t131 = t124 * t113 - t129 * t126;
	t132 = sin(pkin(16));
	t134 = t122 * t121 / 0.4e1 + t132 * t131 / 0.4e1;
	t135 = t134 * t109;
	t136 = t134 ^ 2;
	t139 = -t122 * t131 / 0.4e1 + t132 * t121 / 0.4e1;
	t140 = t139 ^ 2;
	t141 = t136 + t140;
	t142 = sqrt(t141);
	t143 = 0.1e1 / t142;
	t147 = t100 * t104 - t107 * t70;
	t148 = t139 * t147;
	t150 = t143 * t135 - t143 * t148;
	t151 = sin(qJ(3));
	t153 = sin(qJ(1));
	t154 = cos(qJ(3));
	t156 = t151 * t150 - t154 * t153;
	t157 = t61 * t5;
	t159 = -t68 * t2 + t157;
	t162 = t61 * t2;
	t163 = -t68 * t5 - t162;
	t165 = t100 * t159 + t107 * t163;
	t166 = t134 * t165;
	t170 = t100 * t163 - t107 * t159;
	t171 = t139 * t170;
	t173 = t143 * t166 - t143 * t171;
	t174 = 0.1e1 / t173;
	t175 = t174 * t156;
	t176 = 0.1e1 / t151;
	t177 = t2 * t153;
	t178 = t61 * t177;
	t179 = t5 * t153;
	t181 = -t68 * t179 - t178;
	t184 = t61 * t179;
	t185 = t68 * t177 - t184;
	t187 = t100 * t181 + t107 * t185;
	t188 = t134 * t187;
	t192 = t100 * t185 - t107 * t181;
	t193 = t139 * t192;
	t195 = t143 * t188 - t143 * t193;
	t197 = t154 * t1;
	t198 = t151 * t195 + t197;
	t199 = t198 ^ 2;
	t200 = t173 ^ 2;
	t201 = 0.1e1 / t200;
	t203 = t151 ^ 2;
	t204 = 0.1e1 / t203;
	t207 = 0.1e1 / (t204 * t201 * t199 + 0.1e1);
	t208 = t207 * t176;
	t210 = -t13 * pkin(2);
	t212 = 0.1e1 / t27;
	t213 = pkin(1) * t9;
	t217 = pkin(1) * pkin(2);
	t219 = t25 * pkin(2) * t213 + t217 * t9 * t21;
	t226 = -t219 * t212 * t10 - 0.2e1 * t217 * t9 * t30 - t27 * t210 - t43;
	t229 = t49 * t38;
	t232 = t217 * t9 * t56;
	t244 = t9 ^ 2;
	t248 = -0.2e1 * pkin(1) * t244 * t31 + t219 * t212 * t30 + t34 * t210 - t28;
	t253 = 0.1e1 / t48 / t39;
	t255 = pkin(2) * t213;
	t266 = (0.2e1 * t226 * t49 * t46 * t36 + 0.2e1 * t248 * t49 * t46 * t44 + 0.4e1 * t255 * t253 * t47 + 0.4e1 * t255 * t253 * t52) / t55 / t54 * t40;
	t277 = -t57 * t226 * t4 - 0.2e1 * t232 * t229 * t37 + t266 * t38 * t37 / 0.2e1 + t57 * t248 * t7 + 0.2e1 * t232 * t229 * t59 - t266 * t38 * t59 / 0.2e1;
	t296 = -t57 * t226 * t7 - 0.2e1 * t232 * t229 * t64 + t266 * t38 * t64 / 0.2e1 - t57 * t248 * t4 - 0.2e1 * t232 * t229 * t66 + t266 * t38 * t66 / 0.2e1;
	t298 = -t277 * t177 + t68 * t177 - t296 * t179 - t184;
	t303 = t296 * t177 - t277 * t179 + t68 * t179 + t178;
	t313 = t143 * t134 * (t100 * t298 + t107 * t303) - t143 * t139 * (t100 * t303 - t107 * t298);
	t319 = -t296 * t2 + t277 * t5 - t68 * t5 - t162;
	t324 = -t277 * t2 + t68 * t2 - t296 * t5 - t157;
	t334 = t143 * t134 * (t100 * t319 + t107 * t324) - t143 * t139 * (t100 * t324 - t107 * t319);
	t337 = t207 * t201 * t198;
	t339 = -t207 * t174 * t313 + t337 * t176 * t334;
	t341 = t151 * t1;
	t342 = t154 * t195 - t341;
	t349 = t207 * t204 * t198 * t154 * t174 - t208 * t174 * t342;
	t350 = t112 * t72;
	t352 = 0.1e1 / t83;
	t358 = -t80 * t79 * t76 + t80 * t127 - t79 * t78 + t80 * t78;
	t366 = t91 ^ 2;
	t367 = 1 / t366;
	t371 = 0.1e1 / (-t367 * t81 * t78 + 0.1e1);
	t382 = t371 / t91 * t73 * pkin(7) * (-t83 * t350 + t358 * t352 * t75 / 0.2e1) - t371 * t367 * t83 * t73 * pkin(7) * (-t91 * t350 - 0.2e1 * t73 * t75);
	t383 = t382 * t71;
	t385 = t382 * t97;
	t387 = t98 * t383 - t95 * t385;
	t391 = t95 * t383 + t98 * t385;
	t398 = 0.1e1 / t111 / t73 * t110;
	t411 = 0.2e1 * t72 * t73;
	t420 = t115 * t398 / 0.2e1 - t358 * t91 * t72 * t352 * t113 / 0.8e1 + t73 * t114 * t113 / 0.2e1 - t119 * t398 / 0.2e1 + t83 * t411 * t113 / 0.4e1 + t358 * t352 * t72 * t117 * t113 / 0.8e1;
	t434 = t72 * t81;
	t446 = -t124 * t398 / 0.2e1 + t91 * t411 * t113 / 0.4e1 - t73 * t118 * t113 / 0.2e1 + t129 * t76 * t398 / 0.2e1 + t434 * t77 * t113 / 0.4e1 - t434 * t126 / 0.4e1 + t72 * t80 * t77 * t126 / 0.4e1 - t72 * t127 * t126 / 0.4e1;
	t448 = t122 * t420 + t132 * t446;
	t456 = -t122 * t446 + t132 * t420;
	t459 = 0.2e1 * (t448 * t134 + t456 * t139) / t142 / t141;
	t471 = t143 * t134 * (t387 * t181 + t391 * t185) + t143 * t448 * t187 - t459 * t188 / 0.2e1 - t143 * t139 * (-t391 * t181 + t387 * t185) - t143 * t456 * t192 + t459 * t193 / 0.2e1;
	t492 = t143 * t134 * (t387 * t159 + t391 * t163) + t143 * t448 * t165 - t459 * t166 / 0.2e1 - t143 * t139 * (-t391 * t159 + t387 * t163) - t143 * t456 * t170 + t459 * t171 / 0.2e1;
	t495 = -t207 * t174 * t471 + t337 * t176 * t492;
	t498 = -t100 * t181 - t107 * t185;
	t503 = -t100 * t185 + t107 * t181;
	t506 = t143 * t134 * t498 - t143 * t139 * t503;
	t509 = t151 * t173;
	t510 = atan2(t198, -t509);
	t511 = cos(t510);
	t512 = t173 * t511;
	t514 = sin(t510);
	t516 = -t151 * t512 + t198 * t514;
	t517 = 0.1e1 / t516;
	t519 = t156 ^ 2;
	t520 = t516 ^ 2;
	t521 = 0.1e1 / t520;
	t524 = 0.1e1 / (t521 * t519 + 0.1e1);
	t535 = t524 * t521;
	t541 = -t277 * t3 - t296 * t63 + t68 * t3 - t103;
	t546 = -t277 * t63 + t296 * t3 + t68 * t63 + t62;
	t548 = t100 * t541 + t107 * t546;
	t553 = t100 * t546 - t107 * t541;
	t556 = t143 * t134 * t548 - t143 * t139 * t553;
	t558 = t524 * t517;
	t574 = -t154 * t150 - t151 * t153;
	t589 = t391 * t104 + t387 * t70;
	t598 = t387 * t104 - t391 * t70;
	t605 = t143 * t134 * t589 + t143 * t448 * t109 - t459 * t135 / 0.2e1 - t143 * t139 * t598 - t143 * t456 * t147 + t459 * t148 / 0.2e1;
	t621 = t154 * t506 + t341;
	t622 = cos(qJ(4));
	t628 = -t143 * t134 * t503 - t143 * t139 * t498;
	t629 = sin(qJ(4));
	t633 = t139 * t109;
	t635 = t134 * t147;
	t637 = -t143 * t633 - t143 * t635;
	t639 = -t629 * t574 + t622 * t637;
	t640 = 0.1e1 / t639;
	t644 = t622 * t574 + t629 * t637;
	t645 = t644 ^ 2;
	t646 = t639 ^ 2;
	t647 = 0.1e1 / t646;
	t650 = 0.1e1 / (t647 * t645 + 0.1e1);
	t656 = t650 * t647;
	t659 = t154 * t556;
	t665 = -t143 * t134 * t553 - t143 * t139 * t548;
	t687 = t154 * t605;
	t701 = -t143 * t139 * t589 - t143 * t456 * t109 + t459 * t633 / 0.2e1 - t143 * t134 * t598 - t143 * t448 * t147 + t459 * t635 / 0.2e1;
	unknown(1,1) = -t208 * t175;
	unknown(1,2) = t339;
	unknown(1,3) = t349;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t495;
	unknown(2,1) = t524 * t517 * (-t151 * t506 + t197) + t535 * t156 * (-t198 * t511 * t207 * t176 * t175 - t514 * t207 * t156 + t156 * t514);
	unknown(2,2) = -t558 * t151 * t556 + t535 * t156 * (t151 * t313 * t514 - t151 * t334 * t511 + t198 * t511 * t339 + t509 * t514 * t339);
	unknown(2,3) = t524 * t517 * t574 + t535 * t156 * (t198 * t511 * t349 + t509 * t514 * t349 - t154 * t512 + t342 * t514);
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = -t558 * t151 * t605 + t535 * t156 * (t151 * t471 * t514 - t151 * t492 * t511 + t198 * t511 * t495 + t509 * t514 * t495);
	unknown(3,1) = t650 * t640 * (-t622 * t621 + t629 * t628) - t656 * t644 * (t629 * t621 + t622 * t628);
	unknown(3,2) = t650 * t640 * (-t622 * t659 + t629 * t665) - t656 * t644 * (t622 * t665 + t629 * t659);
	unknown(3,3) = t650 * t647 * t644 * t629 * t156 + t650 * t640 * t622 * t156;
	unknown(3,4) = t656 * t644 ^ 2 + t650;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = t650 * t640 * (-t622 * t687 + t629 * t701) - t656 * t644 * (t622 * t701 + t629 * t687);
	Ja_rot = unknown;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:48
	% EndTime: 2020-04-11 12:37:58
	% DurationCPUTime: 32.82s
	% Computational Cost: add. (1794070->255), mult. (2379587->623), div. (183641->32), fcn. (1342158->30), ass. (0->271)
	unknown=NaN(3,6);
	t1 = cos(qJ(1));
	t2 = sin(qJ(2));
	t3 = t2 * t1;
	t4 = cos(pkin(15));
	t5 = cos(qJ(2));
	t7 = sin(pkin(15));
	t9 = -t7 * t2 + t4 * t5;
	t10 = t9 * pkin(2);
	t13 = t4 * t2 + t7 * t5;
	t16 = 0.2e1 * pkin(2) * pkin(1) * t13;
	t17 = pkin(1) ^ 2;
	t21 = -t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4));
	t25 = -t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4));
	t27 = sqrt(-t25 * t21);
	t28 = t27 * t10;
	t30 = -t13 * pkin(2) + pkin(1);
	t31 = pkin(2) ^ 2;
	t32 = pkin(5) ^ 2;
	t33 = pkin(4) ^ 2;
	t34 = -t16 + t17 + t31 - t32 + t33;
	t36 = t34 * t30 - t28;
	t37 = t36 * t4;
	t38 = 0.1e1 / pkin(4);
	t39 = -t16 + t17 + t31;
	t40 = 0.1e1 / t39;
	t43 = t34 * t10;
	t44 = t27 * t30 + t43;
	t45 = t44 ^ 2;
	t46 = 0.1e1 / t33;
	t47 = t46 * t45;
	t48 = t39 ^ 2;
	t49 = 0.1e1 / t48;
	t51 = t36 ^ 2;
	t52 = t46 * t51;
	t54 = t49 * t47 + t49 * t52;
	t55 = sqrt(t54);
	t56 = 0.1e1 / t55;
	t57 = t56 * t40 * t38;
	t59 = t44 * t7;
	t61 = -t57 * t37 + t57 * t59;
	t62 = t61 * t3;
	t63 = t5 * t1;
	t64 = t36 * t7;
	t66 = t44 * t4;
	t68 = -t57 * t64 - t57 * t66;
	t70 = -t68 * t63 - t62;
	t71 = cos(pkin(17));
	t72 = 0.1e1 / pkin(7);
	t73 = qJ(6) + pkin(8);
	t75 = 0.1e1 / t73 * t72;
	t76 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
	t77 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
	t78 = t77 * t76;
	t79 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
	t80 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
	t81 = t80 * t79;
	t83 = sqrt(-t81 * t78);
	t85 = (pkin(6) ^ 2);
	t86 = (pkin(7) ^ 2);
	t87 = (pkin(8) ^ 2);
	t89 = 2 * pkin(8) * qJ(6);
	t90 = qJ(6) ^ 2;
	t91 = t85 - t86 - t87 - t89 - t90;
	t93 = atan2(t83 * t75, t91 * t75);
	t94 = t93 + pkin(16);
	t95 = cos(t94);
	t97 = sin(pkin(17));
	t98 = sin(t94);
	t100 = -t95 * t71 - t98 * t97;
	t103 = t61 * t63;
	t104 = t68 * t3 - t103;
	t107 = t98 * t71 - t95 * t97;
	t109 = t100 * t70 + t107 * t104;
	t110 = 0.1e1 / pkin(6);
	t111 = t73 ^ 2;
	t112 = 0.1e1 / t111;
	t113 = t112 * t110;
	t114 = t72 * t83;
	t115 = t91 * t114;
	t117 = t85 - t86 + t87 + t89 + t90;
	t118 = t72 * t117;
	t119 = t83 * t118;
	t121 = -t115 * t113 + t119 * t113;
	t122 = cos(pkin(16));
	t124 = t91 * t118;
	t126 = t76 * t113;
	t127 = t79 * t77;
	t129 = t72 * t80 * t127;
	t131 = t124 * t113 - t129 * t126;
	t132 = sin(pkin(16));
	t134 = t122 * t121 / 0.4e1 + t132 * t131 / 0.4e1;
	t135 = t134 * t109;
	t136 = t134 ^ 2;
	t139 = -t122 * t131 / 0.4e1 + t132 * t121 / 0.4e1;
	t140 = t139 ^ 2;
	t141 = t136 + t140;
	t142 = sqrt(t141);
	t143 = 0.1e1 / t142;
	t147 = t100 * t104 - t107 * t70;
	t148 = t139 * t147;
	t150 = t143 * t135 - t143 * t148;
	t151 = cos(qJ(3));
	t153 = sin(qJ(1));
	t154 = sin(qJ(3));
	t156 = t151 * t150 + t154 * t153;
	t157 = sin(qJ(4));
	t159 = t139 * t109;
	t161 = t134 * t147;
	t163 = -t143 * t159 - t143 * t161;
	t164 = cos(qJ(4));
	t166 = -t157 * t156 - t164 * t163;
	t167 = t61 * t5;
	t169 = -t68 * t2 + t167;
	t172 = t61 * t2;
	t173 = -t68 * t5 - t172;
	t175 = t100 * t169 + t107 * t173;
	t176 = t134 * t175;
	t180 = t100 * t173 - t107 * t169;
	t181 = t139 * t180;
	t183 = t143 * t176 - t143 * t181;
	t184 = t151 * t183;
	t186 = t139 * t175;
	t188 = t134 * t180;
	t190 = -t143 * t186 - t143 * t188;
	t192 = t157 * t184 + t164 * t190;
	t193 = 0.1e1 / t192;
	t194 = t193 * t166;
	t195 = t2 * t153;
	t196 = t61 * t195;
	t197 = t5 * t153;
	t199 = -t68 * t197 - t196;
	t202 = t61 * t197;
	t203 = t68 * t195 - t202;
	t205 = t100 * t199 + t107 * t203;
	t206 = t134 * t205;
	t210 = t100 * t203 - t107 * t199;
	t211 = t139 * t210;
	t213 = t143 * t206 - t143 * t211;
	t215 = t154 * t1;
	t216 = t151 * t213 - t215;
	t218 = t139 * t205;
	t220 = t134 * t210;
	t222 = -t143 * t218 - t143 * t220;
	t224 = -t157 * t216 - t164 * t222;
	t225 = t224 ^ 2;
	t226 = t192 ^ 2;
	t227 = 0.1e1 / t226;
	t230 = 0.1e1 / (t227 * t225 + 0.1e1);
	t232 = -t13 * pkin(2);
	t234 = 0.1e1 / t27;
	t235 = pkin(1) * t9;
	t239 = pkin(1) * pkin(2);
	t241 = t25 * pkin(2) * t235 + t239 * t9 * t21;
	t248 = -t241 * t234 * t10 - 0.2e1 * t239 * t9 * t30 - t27 * t232 - t43;
	t251 = t49 * t38;
	t254 = t239 * t9 * t56;
	t266 = t9 ^ 2;
	t270 = -0.2e1 * pkin(1) * t266 * t31 + t241 * t234 * t30 + t34 * t232 - t28;
	t275 = 0.1e1 / t48 / t39;
	t277 = pkin(2) * t235;
	t288 = (0.2e1 * t248 * t49 * t46 * t36 + 0.2e1 * t270 * t49 * t46 * t44 + 0.4e1 * t277 * t275 * t47 + 0.4e1 * t277 * t275 * t52) / t55 / t54 * t40;
	t299 = -t57 * t248 * t4 - 0.2e1 * t254 * t251 * t37 + t288 * t38 * t37 / 0.2e1 + t57 * t270 * t7 + 0.2e1 * t254 * t251 * t59 - t288 * t38 * t59 / 0.2e1;
	t318 = -t57 * t248 * t7 - 0.2e1 * t254 * t251 * t64 + t288 * t38 * t64 / 0.2e1 - t57 * t270 * t4 - 0.2e1 * t254 * t251 * t66 + t288 * t38 * t66 / 0.2e1;
	t320 = -t299 * t195 + t68 * t195 - t318 * t197 - t202;
	t325 = t318 * t195 - t299 * t197 + t68 * t197 + t196;
	t327 = t100 * t320 + t107 * t325;
	t332 = t100 * t325 - t107 * t320;
	t344 = -t157 * t151 * (t143 * t134 * t327 - t143 * t139 * t332) - t164 * (-t143 * t134 * t332 - t143 * t139 * t327);
	t350 = -t318 * t2 + t299 * t5 - t68 * t5 - t172;
	t355 = -t299 * t2 + t68 * t2 - t318 * t5 - t167;
	t357 = t100 * t350 + t107 * t355;
	t362 = t100 * t355 - t107 * t350;
	t374 = t157 * t151 * (t143 * t134 * t357 - t143 * t139 * t362) + t164 * (-t143 * t134 * t362 - t143 * t139 * t357);
	t376 = t230 * t227;
	t378 = t230 * t193 * t344 - t376 * t224 * t374;
	t380 = t151 * t1;
	t381 = -t154 * t213 - t380;
	t390 = t230 * t227 * t224 * t157 * t154 * t183 - t230 * t193 * t157 * t381;
	t393 = t157 * t222 - t164 * t216;
	t398 = -t157 * t190 + t164 * t184;
	t401 = t230 * t193 * t393 - t376 * t224 * t398;
	t402 = t112 * t72;
	t404 = 0.1e1 / t83;
	t410 = -t80 * t79 * t76 + t80 * t127 - t79 * t78 + t80 * t78;
	t418 = t91 ^ 2;
	t419 = 1 / t418;
	t423 = 0.1e1 / (-t419 * t81 * t78 + 0.1e1);
	t434 = t423 / t91 * t73 * pkin(7) * (-t83 * t402 + t410 * t404 * t75 / 0.2e1) - t423 * t419 * t83 * t73 * pkin(7) * (-t91 * t402 - 0.2e1 * t73 * t75);
	t435 = t434 * t71;
	t437 = t434 * t97;
	t439 = t98 * t435 - t95 * t437;
	t443 = t95 * t435 + t98 * t437;
	t445 = t439 * t199 + t443 * t203;
	t450 = 0.1e1 / t111 / t73 * t110;
	t463 = 0.2e1 * t72 * t73;
	t472 = t115 * t450 / 0.2e1 - t410 * t91 * t72 * t404 * t113 / 0.8e1 + t73 * t114 * t113 / 0.2e1 - t119 * t450 / 0.2e1 + t83 * t463 * t113 / 0.4e1 + t410 * t404 * t72 * t117 * t113 / 0.8e1;
	t486 = t72 * t81;
	t498 = -t124 * t450 / 0.2e1 + t91 * t463 * t113 / 0.4e1 - t73 * t118 * t113 / 0.2e1 + t129 * t76 * t450 / 0.2e1 + t486 * t77 * t113 / 0.4e1 - t486 * t126 / 0.4e1 + t72 * t80 * t77 * t126 / 0.4e1 - t72 * t127 * t126 / 0.4e1;
	t500 = t122 * t472 + t132 * t498;
	t508 = -t122 * t498 + t132 * t472;
	t511 = 0.2e1 * (t500 * t134 + t508 * t139) / t142 / t141;
	t516 = -t443 * t199 + t439 * t203;
	t540 = -t157 * t151 * (t143 * t134 * t445 + t143 * t500 * t205 - t511 * t206 / 0.2e1 - t143 * t139 * t516 - t143 * t508 * t210 + t511 * t211 / 0.2e1) - t164 * (-t143 * t139 * t445 - t143 * t508 * t205 + t511 * t218 / 0.2e1 - t143 * t134 * t516 - t143 * t500 * t210 + t511 * t220 / 0.2e1);
	t545 = t439 * t169 + t443 * t173;
	t554 = -t443 * t169 + t439 * t173;
	t578 = t157 * t151 * (t143 * t134 * t545 + t143 * t500 * t175 - t511 * t176 / 0.2e1 - t143 * t139 * t554 - t143 * t508 * t180 + t511 * t181 / 0.2e1) + t164 * (-t143 * t139 * t545 - t143 * t508 * t175 + t511 * t186 / 0.2e1 - t143 * t134 * t554 - t143 * t500 * t180 + t511 * t188 / 0.2e1);
	t581 = t230 * t193 * t540 - t376 * t224 * t578;
	t584 = -t100 * t199 - t107 * t203;
	t589 = -t100 * t203 + t107 * t199;
	t592 = t143 * t134 * t584 - t143 * t139 * t589;
	t594 = t151 * t592 + t215;
	t600 = -t143 * t134 * t589 - t143 * t139 * t584;
	t603 = atan2(t224, t192);
	t604 = cos(t603);
	t606 = sin(t603);
	t608 = t192 * t604 + t224 * t606;
	t609 = 0.1e1 / t608;
	t611 = t166 ^ 2;
	t612 = t608 ^ 2;
	t613 = 0.1e1 / t612;
	t616 = 0.1e1 / (t613 * t611 + 0.1e1);
	t626 = t616 * t613;
	t632 = -t299 * t3 + t68 * t3 - t318 * t63 - t103;
	t637 = -t299 * t63 + t318 * t3 + t68 * t63 + t62;
	t639 = t100 * t632 + t107 * t637;
	t644 = t100 * t637 - t107 * t632;
	t647 = t143 * t134 * t639 - t143 * t139 * t644;
	t648 = t151 * t647;
	t654 = -t143 * t134 * t644 - t143 * t139 * t639;
	t671 = -t154 * t150 + t151 * t153;
	t690 = t164 * t156 - t157 * t163;
	t705 = t443 * t104 + t439 * t70;
	t714 = t439 * t104 - t443 * t70;
	t721 = t143 * t134 * t705 + t143 * t500 * t109 - t511 * t135 / 0.2e1 - t143 * t139 * t714 - t143 * t508 * t147 + t511 * t148 / 0.2e1;
	t722 = t151 * t721;
	t736 = -t143 * t139 * t705 - t143 * t508 * t109 + t511 * t159 / 0.2e1 - t143 * t134 * t714 - t143 * t500 * t147 + t511 * t161 / 0.2e1;
	t753 = -t157 * t600 + t164 * t594;
	t754 = sin(qJ(5));
	t757 = -t154 * t592 + t380;
	t758 = cos(qJ(5));
	t763 = t754 * t671 + t758 * t690;
	t764 = 0.1e1 / t763;
	t768 = -t758 * t671 + t754 * t690;
	t769 = t768 ^ 2;
	t770 = t763 ^ 2;
	t771 = 0.1e1 / t770;
	t774 = 0.1e1 / (t771 * t769 + 0.1e1);
	t780 = t774 * t771;
	t785 = -t157 * t654 + t164 * t648;
	t787 = t154 * t647;
	t798 = t164 * t671;
	t823 = -t157 * t736 + t164 * t722;
	t825 = t154 * t721;
	unknown(1,1) = t230 * t194;
	unknown(1,2) = t378;
	unknown(1,3) = t390;
	unknown(1,4) = t401;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = t581;
	unknown(2,1) = t616 * t609 * (t157 * t594 + t164 * t600) + t626 * t166 * (t224 * t604 * t230 * t194 - t606 * t230 * t166 + t166 * t606);
	unknown(2,2) = t616 * t609 * (t157 * t648 + t164 * t654) + t626 * t166 * (-t192 * t606 * t378 + t224 * t604 * t378 + t344 * t606 + t374 * t604);
	unknown(2,3) = t616 * t609 * t157 * t671 + t626 * t166 * (-t157 * t154 * t183 * t604 - t157 * t381 * t606 - t192 * t606 * t390 + t224 * t604 * t390);
	unknown(2,4) = t616 * t609 * t690 + t626 * t166 * (-t192 * t606 * t401 + t224 * t604 * t401 + t393 * t606 + t398 * t604);
	unknown(2,5) = 0.0e0;
	unknown(2,6) = t616 * t609 * (t157 * t722 + t164 * t736) + t626 * t166 * (-t192 * t606 * t581 + t224 * t604 * t581 + t540 * t606 + t578 * t604);
	unknown(3,1) = t774 * t764 * (t754 * t753 - t758 * t757) - t780 * t768 * (t758 * t753 + t754 * t757);
	unknown(3,2) = t774 * t764 * (t754 * t785 + t758 * t787) - t780 * t768 * (-t754 * t787 + t758 * t785);
	unknown(3,3) = t774 * t764 * (t758 * t156 + t754 * t798) - t780 * t768 * (-t754 * t156 + t758 * t798);
	unknown(3,4) = -t774 * t771 * t768 * t758 * t166 + t774 * t764 * t754 * t166;
	unknown(3,5) = t780 * t768 ^ 2 + t774;
	unknown(3,6) = t774 * t764 * (t754 * t823 + t758 * t825) - t780 * t768 * (-t754 * t825 + t758 * t823);
	Ja_rot = unknown;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobia_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:15
	% EndTime: 2020-04-11 12:36:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobia_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:07
	% EndTime: 2020-04-11 12:36:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 11
	%% Symbolic Calculation
	% From jacobia_rot_11_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:28
	% EndTime: 2020-04-11 12:36:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
elseif link_index == 12
	%% Symbolic Calculation
	% From jacobia_rot_12_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 12:36:09
	% EndTime: 2020-04-11 12:36:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->18)
	unknown=NaN(3,6);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	Ja_rot = unknown;
else
	Ja_rot=NaN(3,6);
end