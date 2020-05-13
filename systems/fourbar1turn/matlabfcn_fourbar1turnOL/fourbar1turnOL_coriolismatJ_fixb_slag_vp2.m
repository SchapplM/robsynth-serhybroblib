% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnOL_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:40:50
% DurationCPUTime: 0.20s
% Computational Cost: add. (259->36), mult. (577->56), div. (0->0), fcn. (567->6), ass. (0->23)
t40 = sin(qJ(3));
t41 = sin(qJ(2));
t42 = cos(qJ(3));
t43 = cos(qJ(2));
t31 = t40 * t41 - t42 * t43;
t32 = t40 * t43 + t42 * t41;
t39 = Ifges(4,5) * t31 + Ifges(4,6) * t32;
t47 = t39 * qJD(3);
t28 = t32 ^ 2 * Ifges(4,4) + (-Ifges(4,4) * t31 + (Ifges(4,1) - Ifges(4,2)) * t32) * t31;
t34 = t41 * t43;
t35 = t43 * (-t32 * mrSges(4,1) + t31 * mrSges(4,2));
t1 = (t41 ^ 2 - t43 ^ 2) * Ifges(3,4) + t28 + (-Ifges(3,1) + Ifges(3,2)) * t34 + (-t41 * (-t31 * mrSges(4,1) - t32 * mrSges(4,2)) + t35 + m(4) * pkin(2) * t34) * pkin(2);
t38 = t1 * qJD(1);
t2 = pkin(2) * t35 + t28;
t37 = t2 * qJD(1);
t26 = sin(qJ(4));
t27 = cos(qJ(4));
t4 = (pkin(1) * mrSges(5,2) - Ifges(5,4) * t27) * t27 + (pkin(1) * mrSges(5,1) + Ifges(5,4) * t26 + (-Ifges(5,1) + Ifges(5,2)) * t27) * t26;
t36 = t4 * qJD(1);
t25 = (mrSges(4,1) * t40 + mrSges(4,2) * t42) * pkin(2);
t33 = t25 * qJD(2);
t24 = t25 * qJD(3);
t3 = [-t1 * qJD(2) - t2 * qJD(3) - t4 * qJD(4), -t38 + (Ifges(3,5) * t43 - Ifges(3,6) * t41 + (t42 * t31 - t40 * t32) * mrSges(4,3) * pkin(2) + t39) * qJD(2) + t47, qJD(2) * t39 - t37 + t47, -t36 + (Ifges(5,5) * t27 - Ifges(5,6) * t26) * qJD(4), 0; t38, t24, t24 + t33, 0, 0; t37, -t33, 0, 0, 0; t36, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
Cq = t3;
