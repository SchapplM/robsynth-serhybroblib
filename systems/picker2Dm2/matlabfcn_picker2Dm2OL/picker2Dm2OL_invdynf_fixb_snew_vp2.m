% Calculate vector of cutting forces with Newton-Euler
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% qJDD [12x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x11]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = picker2Dm2OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2OL_invdynf_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:59
% EndTime: 2020-05-09 23:19:03
% DurationCPUTime: 0.76s
% Computational Cost: add. (7942->139), mult. (9941->177), div. (0->0), fcn. (5472->22), ass. (0->91)
t100 = -m(4) - m(10);
t99 = -m(11) - m(5);
t78 = sin(qJ(1));
t87 = cos(qJ(1));
t96 = -t78 * g(1) + t87 * g(2);
t40 = qJDD(1) * pkin(1) + t96;
t90 = qJD(1) ^ 2;
t97 = t87 * g(1) + t78 * g(2);
t41 = -t90 * pkin(1) + t97;
t77 = sin(qJ(2));
t86 = cos(qJ(2));
t31 = t86 * t40 - t77 * t41;
t63 = qJDD(1) + qJDD(2);
t27 = t63 * pkin(3) + t31;
t32 = t77 * t40 + t86 * t41;
t65 = qJD(1) + qJD(2);
t61 = t65 ^ 2;
t29 = -t61 * pkin(3) + t32;
t75 = sin(qJ(4));
t84 = cos(qJ(4));
t98 = t75 * t27 + t84 * t29;
t95 = t84 * t27 - t75 * t29;
t28 = t63 * pkin(2) + t31;
t30 = -t61 * pkin(2) + t32;
t76 = sin(qJ(3));
t85 = cos(qJ(3));
t94 = -t85 * t28 + t76 * t30;
t56 = qJD(3) + t65;
t55 = qJD(4) + t65;
t53 = qJDD(3) + t63;
t52 = qJDD(4) + t63;
t93 = -m(3) - m(7) + t99 + t100;
t92 = -t76 * t28 - t85 * t30;
t91 = -m(2) - m(9) + t93;
t89 = qJD(5) ^ 2;
t88 = qJD(7) ^ 2;
t83 = cos(qJ(5));
t82 = cos(qJ(6));
t81 = cos(qJ(7));
t80 = cos(qJ(8));
t79 = cos(qJ(9));
t74 = sin(qJ(5));
t73 = sin(qJ(6));
t72 = sin(qJ(7));
t71 = sin(qJ(8));
t70 = sin(qJ(9));
t69 = cos(pkin(8));
t68 = cos(qJ(10));
t67 = sin(pkin(8));
t66 = sin(qJ(10));
t64 = qJD(1) + qJD(8);
t62 = qJDD(1) + qJDD(8);
t60 = t64 ^ 2;
t54 = qJD(6) + t65;
t51 = qJDD(6) + t63;
t50 = t56 ^ 2;
t49 = t55 ^ 2;
t48 = t54 ^ 2;
t47 = qJD(9) + t56;
t46 = qJD(10) + t55;
t45 = t47 ^ 2;
t44 = qJDD(9) + t53;
t43 = t46 ^ 2;
t42 = qJDD(10) + t52;
t39 = -t67 * t74 + t69 * t83;
t38 = t67 * t83 + t69 * t74;
t34 = m(8) * (-t81 * g(1) - t72 * g(2)) + qJDD(7) * mrSges(8,1) - t88 * mrSges(8,2);
t33 = m(8) * (-t72 * g(1) + t81 * g(2)) - qJDD(7) * mrSges(8,2) - t88 * mrSges(8,1);
t26 = m(6) * (-t39 * g(1) - t38 * g(2)) - qJDD(5) * mrSges(6,2) - t89 * mrSges(6,1);
t25 = m(6) * (t38 * g(1) - t39 * g(2)) + qJDD(5) * mrSges(6,1) - t89 * mrSges(6,2);
t20 = m(9) * (-t71 * t40 - t80 * t41) - t62 * mrSges(9,2) - t60 * mrSges(9,1);
t19 = m(9) * (-t80 * t40 + t71 * t41) + t62 * mrSges(9,1) - t60 * mrSges(9,2);
t18 = m(7) * (-t73 * t31 - t82 * t32) - t51 * mrSges(7,2) - t48 * mrSges(7,1);
t17 = m(7) * (-t82 * t31 + t73 * t32) + t51 * mrSges(7,1) - t48 * mrSges(7,2);
t16 = -t50 * pkin(6) + t92;
t15 = -t49 * pkin(4) + t98;
t14 = t53 * pkin(6) + t94;
t13 = t52 * pkin(4) + t95;
t12 = m(10) * (-t70 * t14 - t79 * t16) - t44 * mrSges(10,2) - t45 * mrSges(10,1);
t11 = m(10) * (-t79 * t14 + t70 * t16) + t44 * mrSges(10,1) - t45 * mrSges(10,2);
t10 = m(11) * (-t66 * t13 - t68 * t15) - t42 * mrSges(11,2) - t43 * mrSges(11,1);
t9 = m(11) * (-t68 * t13 + t66 * t15) + t42 * mrSges(11,1) - t43 * mrSges(11,2);
t8 = m(4) * t92 - t50 * mrSges(4,1) - t53 * mrSges(4,2) + t70 * t11 - t79 * t12;
t7 = m(4) * t94 + t53 * mrSges(4,1) - t50 * mrSges(4,2) - t79 * t11 - t70 * t12;
t6 = m(5) * t98 - t49 * mrSges(5,1) - t52 * mrSges(5,2) - t68 * t10 + t66 * t9;
t5 = m(5) * t95 + t52 * mrSges(5,1) - t49 * mrSges(5,2) - t66 * t10 - t68 * t9;
t4 = m(3) * t32 - t61 * mrSges(3,1) - t63 * mrSges(3,2) + t73 * t17 - t82 * t18 - t75 * t5 + t84 * t6 + t76 * t7 - t85 * t8;
t3 = m(3) * t31 + t63 * mrSges(3,1) - t61 * mrSges(3,2) - t82 * t17 - t73 * t18 + t84 * t5 + t75 * t6 - t85 * t7 - t76 * t8;
t2 = m(2) * t97 - t90 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t71 * t19 - t80 * t20 - t77 * t3 + t86 * t4;
t1 = m(2) * t96 + qJDD(1) * mrSges(2,1) - t90 * mrSges(2,2) - t80 * t19 - t71 * t20 + t86 * t3 + t77 * t4;
t21 = [-m(1) * g(1) + t78 * t1 - t87 * t2 - t38 * t25 + t39 * t26 + t72 * t33 + t81 * t34, t2, t4, t8, t6, t26, t18, t33, t20, t12, t10, 0, 0, 0, 0, 0; -m(1) * g(2) - t87 * t1 - t78 * t2 + t39 * t25 + t38 * t26 - t81 * t33 + t72 * t34, t1, t3, t7, t5, t25, t17, t34, t19, t11, t9, 0, 0, 0, 0, 0; (-m(1) - m(6) - m(8) + t91) * g(3), t91 * g(3), t93 * g(3), t100 * g(3), t99 * g(3), -m(6) * g(3), -m(7) * g(3), -m(8) * g(3), -m(9) * g(3), -m(10) * g(3), -m(11) * g(3), 0, 0, 0, 0, 0;];
f_new = t21;
