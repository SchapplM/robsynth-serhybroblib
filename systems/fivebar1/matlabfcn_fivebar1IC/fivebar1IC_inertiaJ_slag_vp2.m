% Calculate joint inertia matrix for
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fivebar1IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1IC_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1IC_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:11
% EndTime: 2020-04-27 06:19:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (172->38), mult. (234->57), div. (34->7), fcn. (96->11), ass. (0->30)
t55 = sin(qJ(2));
t77 = pkin(2) * t55;
t54 = sin(qJ(4));
t76 = pkin(3) * t54;
t75 = Ifges(3,3) / pkin(5) ^ 2;
t74 = Ifges(5,3) / pkin(4) ^ 2;
t69 = qJ(4) + qJ(3);
t70 = qJ(1) + qJ(2);
t71 = cos(0.2e1 * t69) - cos(0.2e1 * t70);
t44 = 0.1e1 / t71;
t58 = 0.2e1 * qJ(3);
t59 = 0.2e1 * qJ(1);
t73 = (t71 * pkin(4) + (cos(qJ(4) + t58) - cos(qJ(4) + t59 + 0.2e1 * qJ(2))) * pkin(3)) * t44;
t72 = (t71 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t58) + cos(qJ(2) + t59)) * pkin(2)) * t44;
t68 = mrSges(3,1) * pkin(2) * cos(qJ(2));
t62 = 0.1e1 / pkin(4);
t67 = t62 * t73;
t60 = 0.1e1 / pkin(5);
t66 = t60 * t72;
t52 = mrSges(3,2) * t77;
t46 = Ifges(3,3) + t52 - t68;
t53 = mrSges(5,1) * pkin(3) * cos(qJ(4));
t45 = -mrSges(5,2) * t76 + Ifges(5,3) + t53;
t64 = pkin(3) ^ 2;
t51 = -sin(t69 - t70);
t48 = 0.1e1 / t51 ^ 2;
t47 = 0.1e1 / t51;
t41 = -Ifges(3,3) * t66 + t46;
t40 = -Ifges(5,3) * t67 + t45;
t1 = [-0.2e1 * t68 + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t52 + (t48 * t55 ^ 2 * t74 + m(3)) * pkin(2) ^ 2 - (t41 + t46) * t66, (t41 * t60 * t76 + (t45 * t62 - t73 * t74) * t77) * t47; (t40 * t62 * t77 + (t46 * t60 - t72 * t75) * t76) * t47, t64 * m(5) + Ifges(4,3) + Ifges(5,3) + 0.2e1 * t53 + (t48 * t64 * t54 * t75 - 0.2e1 * pkin(3) * mrSges(5,2)) * t54 - (t40 + t45) * t67;];
Mq = t1;
