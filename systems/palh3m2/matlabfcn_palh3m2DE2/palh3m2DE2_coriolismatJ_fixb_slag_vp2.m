% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m2DE2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:22:22
% EndTime: 2020-05-07 04:22:23
% DurationCPUTime: 1.59s
% Computational Cost: add. (2434->234), mult. (2480->296), div. (0->0), fcn. (1310->80), ass. (0->152)
t954 = qJ(2) + qJ(3);
t852 = sin(t954);
t860 = pkin(17) + pkin(18);
t929 = pkin(15) + t860;
t918 = (pkin(16) + t929);
t908 = (qJ(4) + t918);
t902 = qJ(2) + t908;
t899 = qJ(3) + t902;
t909 = (-qJ(4) + t918);
t903 = -qJ(2) + t909;
t900 = -qJ(3) + t903;
t883 = sin(t899) / 0.4e1 - sin(t900) / 0.4e1;
t884 = cos(t899) / 0.4e1 + cos(t900) / 0.4e1;
t921 = 0.2e1 * t954;
t895 = sin(t921) / 0.2e1;
t911 = -qJ(2) + t918;
t904 = -qJ(3) + t911;
t910 = qJ(2) + t918;
t905 = qJ(3) + t910;
t920 = -qJ(2) + t929;
t913 = -qJ(3) + t920;
t919 = qJ(2) + t929;
t914 = qJ(3) + t919;
t976 = m(5) + m(6);
t941 = t976 * pkin(4);
t917 = t941 + mrSges(4,1) + mrSges(9,1);
t924 = -pkin(10) * m(6) + mrSges(5,2) - mrSges(6,3);
t931 = cos(t954);
t937 = pkin(8) * m(6) + mrSges(5,1);
t964 = mrSges(4,2) + mrSges(9,2);
t830 = -qJ(2) + t908;
t810 = -qJ(3) + t830;
t970 = cos(t810);
t831 = qJ(2) + t909;
t809 = qJ(3) + t831;
t971 = cos(t809);
t982 = -sin(t809) / 0.4e1 + sin(t810) / 0.4e1;
t992 = (-Ifges(4,1) - Ifges(9,1) + Ifges(4,2) + Ifges(9,2)) * t895 + ((cos(t913) / 0.2e1 + cos(t914) / 0.2e1) * mrSges(9,2) + (-sin(t913) / 0.2e1 + sin(t914) / 0.2e1) * mrSges(9,1)) * pkin(3) - (Ifges(4,4) + Ifges(9,4)) * cos(t921) + (-t917 * t852 - t964 * t931) * pkin(12) + (t895 * t941 + (-t971 / 0.4e1 - t970 / 0.4e1 + t884) * mrSges(6,2) + (t883 - t982) * mrSges(6,1) + (sin(t905) / 0.2e1 - sin(t904) / 0.2e1) * t937 + (-cos(t905) / 0.2e1 + cos(t904) / 0.2e1) * t924) * pkin(4);
t979 = 2 * qJ(2);
t991 = sin(t979) / 0.2e1;
t868 = sin(qJ(2));
t871 = cos(qJ(3));
t936 = -pkin(4) * t871 + pkin(1);
t867 = sin(qJ(3));
t872 = cos(qJ(2));
t956 = t867 * t872;
t758 = pkin(4) * t956 - t936 * t868;
t957 = t867 * t868;
t838 = pkin(4) * t957;
t759 = t936 * t872 + t838;
t869 = sin(pkin(15));
t873 = cos(pkin(15));
t734 = t758 * t873 - t759 * t869;
t863 = sin(pkin(16));
t864 = cos(pkin(16));
t926 = t758 * t869 + t759 * t873;
t990 = t734 * t863 + t864 * t926;
t989 = t734 * t864 - t863 * t926;
t955 = t871 * t872;
t793 = t955 - t957;
t792 = t868 * t871 + t956;
t958 = t792 * t873;
t738 = t793 * t869 + t958;
t959 = t792 * t869;
t925 = t793 * t873 - t959;
t988 = -t738 * t863 + t864 * t925;
t987 = t738 * t864 + t863 * t925;
t978 = 2 * qJ(4);
t859 = -qJ(2) - pkin(15) + pkin(14);
t977 = 0.2e1 * t859;
t975 = pkin(4) * mrSges(6,1);
t974 = pkin(4) * mrSges(6,2);
t973 = pkin(8) * mrSges(6,1);
t972 = pkin(8) * mrSges(6,2);
t969 = cos(t859);
t966 = pkin(8) * t869;
t965 = pkin(8) * t873;
t866 = sin(qJ(4));
t962 = mrSges(6,1) * t866;
t870 = cos(qJ(4));
t961 = mrSges(6,2) * t870;
t824 = mrSges(6,1) * t870 - mrSges(6,2) * t866;
t960 = t758 * t824;
t862 = qJ(2) + qJ(4);
t953 = qJ(3) + t979;
t952 = qJ(4) - qJ(2);
t951 = qJD(4) * t824;
t791 = t917 * pkin(1) * sin(t953);
t829 = pkin(1) * t964 * cos(t953);
t844 = sin(t859);
t885 = sin(t903) / 0.4e1 - sin(t902) / 0.4e1;
t886 = -cos(t903) / 0.4e1 - cos(t902) / 0.4e1;
t943 = pkin(15) + pkin(18);
t932 = qJ(2) + t943;
t933 = -qJ(2) + t943;
t934 = -cos(t830) / 0.4e1 - cos(t831) / 0.4e1;
t935 = -sin(t830) / 0.4e1 + sin(t831) / 0.4e1;
t721 = -t791 + (-Ifges(3,1) + Ifges(3,2)) * t991 - (Ifges(7,2) - Ifges(7,1)) * sin(t977) / 0.2e1 - Ifges(3,4) * cos(t979) - t829 - Ifges(7,4) * cos(t977) + (mrSges(3,1) * t868 + mrSges(3,2) * t872) * pkin(12) + ((cos(t932) / 0.2e1 - cos(t933) / 0.2e1) * mrSges(8,2) + (sin(t933) / 0.2e1 - sin(t932) / 0.2e1) * mrSges(8,1) + (-sin(t919) / 0.2e1 + sin(t920) / 0.2e1) * m(9) * pkin(3) + (t886 - t934) * mrSges(6,2) + (t885 - t935) * mrSges(6,1) + (-sin(t910) / 0.2e1 + sin(t911) / 0.2e1) * t937 + (pkin(1) * t991 + pkin(12) * t868) * (m(8) + m(9) + m(4) + t976) + (-cos(t911) / 0.2e1 + cos(t910) / 0.2e1) * t924) * pkin(1) + (mrSges(7,1) * t844 - mrSges(7,2) * t969) * pkin(6) + t992;
t950 = t721 * qJD(1);
t855 = pkin(10) * mrSges(6,2) - Ifges(6,6);
t856 = pkin(10) * mrSges(6,1) - Ifges(6,5);
t897 = 2 * t908;
t898 = 2 * t909;
t922 = 2 * t918;
t915 = qJ(4) + t922;
t916 = -qJ(4) + t922;
t722 = -(t855 + t973) * sin(t916) / 0.4e1 - (-t856 - t972) * cos(t915) / 0.4e1 + (t962 / 0.2e1 + t961 / 0.2e1) * pkin(8) + (-t856 + t972) * cos(t916) / 0.4e1 + (-t855 + t973) * sin(t915) / 0.4e1 + (cos(t978) / 0.2e1 - cos(t898) / 0.4e1 - cos(t897) / 0.4e1) * Ifges(6,4) + ((-cos(t908) / 0.2e1 - cos(t909) / 0.2e1) * mrSges(6,2) + (sin(t909) / 0.2e1 - sin(t908) / 0.2e1) * mrSges(6,1)) * pkin(12) + (t883 * mrSges(6,1) + t884 * mrSges(6,2)) * pkin(4) + ((t886 + t934) * mrSges(6,2) + (t885 + t935) * mrSges(6,1)) * pkin(1) + t982 * t975 + (t970 + t971) * t974 / 0.4e1 + (-sin(t898) / 0.8e1 + sin(t897) / 0.8e1 - sin(t978) / 0.4e1) * (Ifges(6,2) - Ifges(6,1));
t949 = t722 * qJD(1);
t906 = t917 * t867;
t930 = t964 * t871;
t725 = -t791 / 0.2e1 - t829 / 0.2e1 + (-t906 / 0.2e1 - t930 / 0.2e1) * pkin(1) + t992;
t948 = t725 * qJD(1);
t850 = -sin(t952);
t851 = sin(t862);
t853 = cos(t952);
t854 = cos(t862);
t857 = qJ(3) + t862;
t842 = sin(t857);
t858 = qJ(3) - t952;
t843 = sin(t858);
t845 = cos(t857);
t846 = cos(t858);
t894 = (-t846 / 0.4e1 + t845 / 0.4e1) * mrSges(6,2) + (t843 / 0.4e1 + t842 / 0.4e1) * mrSges(6,1);
t891 = t894 * pkin(4);
t882 = t891 + ((t853 / 0.4e1 - t854 / 0.4e1) * mrSges(6,2) + (-t850 / 0.4e1 - t851 / 0.4e1) * mrSges(6,1)) * pkin(1);
t729 = t960 / 0.2e1 + t882;
t947 = t729 * qJD(1);
t731 = t852 * (Ifges(4,6) + Ifges(9,6)) + (pkin(4) * mrSges(5,3) - Ifges(4,5) - Ifges(9,5)) * t931 + (t842 / 0.2e1 - t843 / 0.2e1) * t975 + (t845 + t846) * t974 / 0.2e1;
t946 = t731 * qJD(3);
t901 = t792 * t824;
t732 = (t901 / 0.2e1 + t894) * pkin(4);
t945 = t732 * qJD(1);
t823 = t961 + t962;
t848 = sin(t860);
t849 = cos(t860);
t940 = (-t848 * t987 + t988 * t849) * t823 * pkin(4);
t923 = t940 / 0.2e1;
t907 = mrSges(6,1) * t966 + t856 * t873;
t892 = -t873 * t955 + t959;
t893 = t869 * t955 + t958;
t878 = (((-t863 * t869 + t864 * t873) * t849 + t848 * (-t863 * t873 - t864 * t869)) * t838 + ((t893 * t863 + t892 * t864) * t849 + t848 * (-t863 * t892 + t893 * t864)) * pkin(4)) * t823;
t723 = t923 + t878 / 0.2e1;
t750 = (t906 + t930) * pkin(1);
t896 = -qJD(2) * t750 + qJD(4) * t723;
t785 = mrSges(6,2) * t966 + t855 * t873;
t784 = mrSges(6,1) * t965 - t856 * t869;
t783 = mrSges(6,2) * t965 - t855 * t869;
t746 = t750 * qJD(3);
t733 = -pkin(4) * t901 / 0.2e1 + t891;
t730 = -t960 / 0.2e1 + t882;
t724 = -t878 / 0.2e1 + t923;
t1 = [-qJD(2) * t721 - qJD(3) * t725 - qJD(4) * t722, t730 * qJD(4) + t946 - t950 + (Ifges(3,5) * t872 + Ifges(7,5) * t969 - Ifges(3,6) * t868 + Ifges(7,6) * t844 + t731 + ((-mrSges(4,3) - mrSges(5,3) - mrSges(8,3) - mrSges(9,3)) * t872 + (-t853 / 0.2e1 - t854 / 0.2e1) * mrSges(6,2) + (t850 / 0.2e1 - t851 / 0.2e1) * mrSges(6,1)) * pkin(1)) * qJD(2), t731 * qJD(2) + t733 * qJD(4) + t946 - t948, -t949 + t730 * qJD(2) + t733 * qJD(3) + ((-(t784 * t864 - t907 * t863) * t866 + (-t783 * t864 + t785 * t863) * t870) * t849 + (-(-t863 * t784 - t907 * t864) * t866 + (t783 * t863 + t785 * t864) * t870) * t848 + (pkin(12) + t759) * t823) * qJD(4); qJD(4) * t729 + t950, t746, t746 - t896, t947 - t723 * qJD(3) - (t848 * t990 - t989 * t849) * t951; qJD(4) * t732 + t948, t896, 0, t945 + t723 * qJD(2) + (t848 * t988 + t987 * t849) * pkin(4) * t951; -qJD(2) * t729 - qJD(3) * t732 + t949, -t947 - (t848 * t989 + t990 * t849) * t823 * qJD(2) + t724 * qJD(3), qJD(2) * t724 + qJD(3) * t940 - t945, 0;];
Cq = t1;
