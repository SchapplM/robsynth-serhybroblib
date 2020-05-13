% Calculate joint inertia matrix for
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
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
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnIC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnIC_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnIC_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:35
% EndTime: 2020-05-07 11:33:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (105->39), mult. (150->68), div. (16->5), fcn. (137->8), ass. (0->20)
t59 = -qJ(4) + qJ(2);
t47 = sin(qJ(3) + t59);
t46 = 0.1e1 / t47;
t61 = (-pkin(3) * t47 + pkin(2) * sin(t59)) * t46;
t51 = sin(qJ(3));
t52 = sin(qJ(2));
t54 = cos(qJ(3));
t55 = cos(qJ(2));
t44 = t51 * t52 - t54 * t55;
t45 = -t51 * t55 - t54 * t52;
t60 = Ifges(4,5) * t45 + Ifges(4,6) * t44;
t58 = t54 * pkin(2) * mrSges(4,1);
t57 = pkin(2) ^ 2;
t56 = 0.1e1 / pkin(3);
t53 = cos(qJ(4));
t50 = sin(qJ(4));
t49 = t51 ^ 2;
t48 = t51 * pkin(2) * mrSges(4,2);
t40 = Ifges(3,5) * t52 + Ifges(3,6) * t55 + t56 * t60 * t61 + (t51 / pkin(4) * t46 * (Ifges(5,5) * t50 + Ifges(5,6) * t53) + (-t51 * t44 + t54 * t45) * mrSges(4,3)) * pkin(2) + t60;
t1 = [m(5) * pkin(1) ^ 2 + t50 * (Ifges(5,1) * t50 + Ifges(5,4) * t53) + t53 * (Ifges(5,4) * t50 + Ifges(5,2) * t53) - 0.2e1 * pkin(1) * (-t53 * mrSges(5,1) + t50 * mrSges(5,2)) + Ifges(4,1) * t45 ^ 2 + Ifges(3,1) * t52 ^ 2 + Ifges(2,3) + (0.2e1 * Ifges(4,4) * t45 + Ifges(4,2) * t44) * t44 + (-0.2e1 * pkin(2) * (-t44 * mrSges(4,1) + t45 * mrSges(4,2)) + 0.2e1 * Ifges(3,4) * t52 + (m(4) * t57 + Ifges(3,2)) * t55) * t55, t40; t40, -0.2e1 * t58 + Ifges(3,3) + Ifges(4,3) + 0.2e1 * t48 + (m(4) * (t54 ^ 2 + t49) + t49 / pkin(4) ^ 2 / t47 ^ 2 * Ifges(5,3)) * t57 + (0.2e1 * (Ifges(4,3) + t48 - t58) * t56 + Ifges(4,3) * t56 ^ 2 * t61) * t61;];
Mq = t1;
